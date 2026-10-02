/*
    puera - visualizing y-chromosome samples on trees
    Copyright (C) 2025 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "options/tree_table.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/tree/formats/table/reader.hpp"
#include "genesis/tree/function/functions.hpp"
#include "genesis/utils/containers/dataframe.hpp"
#include "genesis/utils/containers/dataframe/operators.hpp"
#include "genesis/utils/containers/dataframe/reader.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>

using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//      Setup Functions
// =================================================================================================

void TreeTableOptions::add_tree_table_opt_to_app(
        CLI::App* sub,
        bool required,
        std::string const& group
) {
    // Correct setup check.
    if( tree_table_opt_.option != nullptr ) {
        throw std::domain_error( "Cannot use the same TreeTableOptions object multiple times." );
    }

    // Input files.
    tree_table_opt_.option = sub->add_option(
        "--tree-table-file",
        tree_table_opt_.value,
        "Provide a table file that describes the tree. The table needs to have (at least) "
        "two columns naming a node and its parent, and contain a row for every node in the tree, "
        "giving its name (in the `--tree-table-id-column`) and the name of its parent (in the "
        "`--tree-table-parent-column`). The root node can have an arbitrary name, and "
        "needs to be the only node that is only ever listed as a parent."
    );
    if( required ) {
        tree_table_opt_.option->required();
    }
    tree_table_opt_.option->check( CLI::ExistingFile );
    tree_table_opt_.option->group( group );

    // Table column names.
    idx_col_opt_.option = sub->add_option(
        "--tree-table-id-column",
        idx_col_opt_.value,
        "Column name for the node id column in the tree table."
    );
    idx_col_opt_.option->group( group );
    par_col_opt_.option = sub->add_option(
        "--tree-table-parent-column",
        par_col_opt_.value,
        "Column name for the parent id column in the tree table."
    );
    par_col_opt_.option->group( group );

    // Separator char
    separator_char_opt_.option = sub->add_option(
        "--tree-table-separator-char",
        separator_char_opt_.value,
        "Separator char between fields of the tree table."
    )->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
    separator_char_opt_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::tree::Tree TreeTableOptions::get_tree() const
{
    read_tree_();
    return tree_;
}

std::unordered_map<std::string, size_t> TreeTableOptions::make_node_name_to_edge_index() const
{
    // Get the tree
    read_tree_();

    // Create a map of branch names to edge indices, for speed
    std::unordered_map<std::string, size_t> result;
    for( auto const& node : tree_.nodes() ) {
        if( is_root( node )) {
            continue;
        }
        auto const& node_name = node.data<CommonNodeData>().name;
        if( result.count( node_name ) > 0 ) {
            throw std::runtime_error( "Duplicate node name: " + node_name );
        }
        auto const edge_index = node.primary_edge().index();
        result[ node_name ] = edge_index;
    }
    return result;
}

void TreeTableOptions::set_tree_branch_length_to_snp_counts(
    std::vector<size_t> const& edge_snp_counts
) const {
    internal_check( edge_snp_counts.size() == tree_.edge_count() );
    for( size_t i = 0; i < edge_snp_counts.size(); ++i ) {
        tree_.edge_at( i ).data<CommonEdgeData>().branch_length = edge_snp_counts[i];
    }
}

// =================================================================================================
//      Internal Functions
// =================================================================================================

// // Custom hash function for std::pair<std::string, std::string>
// struct pair_hash {
//     std::size_t operator()(const std::pair<std::string, std::string>& p) const {
//         auto hash1 = std::hash<std::string>{}(p.first);
//         auto hash2 = std::hash<std::string>{}(p.second);
//         return hash1 ^ (hash2 << 1); // Combine the two hash values
//     }
// };

// // Custom equality function for std::pair<std::string, std::string>
// struct pair_equal {
//     bool operator()(const std::pair<std::string, std::string>& p1,
//                     const std::pair<std::string, std::string>& p2) const {
//         return p1.first == p2.first && p1.second == p2.second;
//     }
// };

void TreeTableOptions::read_tree_() const
{
    if( ! tree_.empty() ) {
        return;
    }

    // LOG_DBG << "read_tree_from_table_file";

    // Read the input table into a dataframe.
    auto const sep_char = translate_separator_char( separator_char_opt_ );
    auto reader = DataframeReader<std::string>( sep_char ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( tree_table_opt_.value ));
    // LOG_DBG << "tree table columns: " << join( table.col_names() );
    auto const& raw_children = table[idx_col_opt_.value].as<std::string>().to_vector();
    auto const& raw_parents  = table[par_col_opt_.value].as<std::string>().to_vector();

    tree_ = make_tree_from_parents_table( raw_children, raw_parents );

    // // Make lists that do not contain duplicates.
    // // Super inefficient, but good enough for now.
    // std::vector<std::string> child_names;
    // std::vector<std::string> parent_names;
    // std::unordered_set<std::pair<std::string, std::string>, pair_hash, pair_equal> duplicates;
    // for( size_t i = 0; i < raw_children.size(); ++i ) {
    //     auto cp = std::make_pair( raw_children[i], raw_parents[i] );
    //     // if( cp.first.empty() || cp.second.empty() ) {
    //     //     continue;
    //     // }
    //     // if( cp.first == "A00" ) {
    //     //     continue;
    //     // }
    //     // if( cp.second.empty() ) {
    //     //     cp.second = "root";
    //     // }
    //     // if( duplicates.count( cp ) > 0 ) {
    //     //     continue;
    //     // }
    //     duplicates.insert( cp );
    //     child_names.push_back( cp.first );
    //     parent_names.push_back( cp.second );
    //
    //     // LOG_DBG << cp.first << " --> " << cp.second;
    // }
    // return make_tree_from_parents_table( child_names, parent_names );
}
