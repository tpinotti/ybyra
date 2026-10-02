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

#include "options/samples.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/tree/format/table/reader.hpp"
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

void SamplesOptions::add_samples_opt_to_app(
        CLI::App* sub,
        bool required,
        std::string const& group
) {
    // Set up the file/dir path
    samples_opt_.add_multi_file_input_opt_to_app(
        sub, "samples", "sample", "calls", "calls", required, group
    );

    // Samples column names.
    idx_col_opt_.option = sub->add_option(
        "--sample-table-id-column",
        idx_col_opt_.value,
        "Column name for the node id column in the sample tables."
    );
    idx_col_opt_.option->group( group );
    par_col_opt_.option = sub->add_option(
        "--sample-table-parent-column",
        par_col_opt_.value,
        "Column name for the parent id column in the sample tables."
    );
    par_col_opt_.option->group( group );
    stt_col_opt_.option = sub->add_option(
        "--sample-table-state-column",
        stt_col_opt_.value,
        "Column name for the state column (ancestral or derived) in the sample tables."
    );
    stt_col_opt_.option->group( group );
    dmg_col_opt_.option = sub->add_option(
        "--sample-table-damage-column",
        dmg_col_opt_.value,
        "Column name for the damage column (yes or no) in the sample tables."
    );
    dmg_col_opt_.option->group( group );

    // Separator char
    separator_char_opt_.option = sub->add_option(
        "--sample-table-separator-char",
        separator_char_opt_.value,
        "Separator char between fields of the sample tables."
    )->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
    separator_char_opt_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::vector<double> SamplesOptions::make_sample_edge_snp_counts_derived(
    TreeTableOptions const& tree_opts,
    size_t smp_idx
) {
    auto const sample_table_infile = samples_opt_.file_path( smp_idx );
    auto const& tree = tree_opts.get_tree();

    // Create a map of branch names to edge indices, for speed
    auto const node_name_to_edge_index = tree_opts.make_node_name_to_edge_index();

    // Read the table
    auto const sep_char = translate_separator_char( separator_char_opt_ );
    auto reader = DataframeReader<std::string>( sep_char );
    reader.row_names_from_first_col( false );
    auto const table = reader.read( from_file( sample_table_infile ));

    // Shortcuts for the columns of the table that we need.
    auto const& col_id     = table[ idx_col_opt_.value ].as<std::string>().to_vector();
    auto const& col_parent = table[ par_col_opt_.value ].as<std::string>().to_vector();
    auto const& col_state  = table[ stt_col_opt_.value ].as<std::string>().to_vector();
    auto const& col_damage = table[ dmg_col_opt_.value ].as<std::string>().to_vector();

    // Prepare the result vector and loop the table
    auto edge_values = std::vector<double>( tree.edge_count(), 0.0 );
    for( size_t i = 0; i < table.rows(); ++i ) {
        if( col_id[i].empty() ) {
            continue;
        }
        if( exclude_damage_.value && col_damage[i] == "yes" ) {
            continue;
        }
        if( col_damage[i] != "yes" && col_damage[i] != "no" ) {
            throw std::runtime_error(
                "Invalid damage column value \"" + col_damage[i] + "\""
            );
        }
        if( node_name_to_edge_index.count( col_id[i] ) == 0 ) {
            throw std::runtime_error(
                "No child with name " + col_id[i]
            );
        }
        auto const edge_index = node_name_to_edge_index.at( col_id[i] );
        auto const& edge = tree.edge_at( edge_index );
        auto const& parent_name = edge.primary_node().data<CommonNodeData>().name;
        auto const& child_name = edge.secondary_node().data<CommonNodeData>().name;
        if( parent_name != col_parent[i] ) {
            throw std::runtime_error(
                "Wrong parent name " + parent_name + " instead of " + col_parent[i]
            );
        }
        if( child_name != col_id[i] ) {
            throw std::runtime_error(
                "Wrong child name " + child_name + " instead of " + col_id[i]
            );
        }

        if( col_state[i] == "derived" || col_state[i] == "+" ) {
            edge_values[edge_index] += 1.0;
        } else if( col_state[i] == "ancestral" || col_state[i] == "-" ) {
            edge_values[edge_index] -= 1.0;
        } else {
            throw std::runtime_error(
                "Invalid state column value \"" + col_state[i] + "\""
            );
        }
    }
    return edge_values;
}
