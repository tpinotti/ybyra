/*
    puera - visualizing y-chromosome samples on trees
    Copyright (C) 2025-2026 Lucas Czech

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

#include "options/snps_table.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/tree/function/function.hpp"
#include "genesis/util/container/dataframe.hpp"
#include "genesis/util/container/dataframe/operator.hpp"
#include "genesis/util/container/dataframe/reader.hpp"
#include "genesis/util/core/fs.hpp"
#include "genesis/util/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>

using namespace genesis::tree;
using namespace genesis::util::container;
using namespace genesis::util::io;

// =================================================================================================
//      Setup Functions
// =================================================================================================

void SnpsTableOptions::add_snps_table_opt_to_app(
        CLI::App* sub,
        bool required,
        std::string const& group
) {
    // Correct setup check.
    if( snps_table_opt_.option != nullptr ) {
        throw std::domain_error( "Cannot use the same SnpsTableOptions object multiple times." );
    }

    // Input files.
    snps_table_opt_.option = sub->add_option(
        "--snps-table-file",
        snps_table_opt_.value,
        "Provide a table file that describes the SNPs. The table needs to have (at least) "
        "two columns naming the node and its parent of the branch where the SNP is located, "
        "and contain a row for every SNP."
    );
    if( required ) {
        snps_table_opt_.option->required();
    }
    snps_table_opt_.option->check( CLI::ExistingFile );
    snps_table_opt_.option->group( group );

    // Table column names.
    idx_col_opt_.option = sub->add_option(
        "--snps-table-id-column",
        idx_col_opt_.value,
        "Column name for the node id column in the SNPs table."
    );
    idx_col_opt_.option->group( group );
    par_col_opt_.option = sub->add_option(
        "--snps-table-parent-column",
        par_col_opt_.value,
        "Column name for the parent id column in the SNPs table."
    );
    par_col_opt_.option->group( group );

    // Separator char
    separator_char_opt_.option = sub->add_option(
        "--snps-table-separator-char",
        separator_char_opt_.value,
        "Separator char between fields of the SNPs table."
    )->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
    separator_char_opt_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::vector<size_t> SnpsTableOptions::get_snps_per_edge_counts(
    TreeTableOptions const& tree_opts
) const {
    // Create a map of branch names to edge indices, for speed
    auto const& tree = tree_opts.get_tree();
    auto const& node_name_to_edge_index = tree_opts.node_name_to_edge_index();

    // First read the table containing the snps and branches
    auto const sep_char = translate_separator_char( separator_char_opt_ );
    auto reader = DataframeReader<std::string>( sep_char ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( snps_table_opt_.value ));

    // Shortcuts for the columns of the table that we need.
    auto const& col_idx = table[idx_col_opt_.value].as<std::string>().to_vector();
    auto const& col_par = table[par_col_opt_.value].as<std::string>().to_vector();
    if( col_idx.size() != col_par.size() ) {
        throw std::runtime_error( "Inconsistent column sizes in SNP table" );
    }

    // Finally, create a vector of edge indices counting the snps on each edge.
    auto edge_snp_counts = std::vector<size_t>( tree.edge_count() );
    for( size_t i = 0; i < col_idx.size(); ++i ) {
        if( col_idx[i].empty() ) {
            LOG_WARN << "Empty SNP table id at position " << i;
            continue;
        }
        if( node_name_to_edge_index.count( col_idx[i] ) == 0 ) {
            throw std::runtime_error(
                "No child with name " + col_idx[i]
            );
        }
        auto const edge_index = node_name_to_edge_index.at( col_idx[i] );
        auto const& edge = tree.edge_at( edge_index );
        auto const& parent_name = edge.primary_node().data<CommonNodeData>().name;
        auto const& child_name  = edge.secondary_node().data<CommonNodeData>().name;
        if( parent_name != col_par[i] ) {
            throw std::runtime_error(
                "Wrong parent name " + parent_name + " instead of " + col_par[i]
            );
        }
        if( child_name != col_idx[i] ) {
            throw std::runtime_error(
                "Wrong child name " + child_name + " instead of " + col_idx[i]
            );
        }
        ++edge_snp_counts[ edge_index ];
    }
    return edge_snp_counts;
}
