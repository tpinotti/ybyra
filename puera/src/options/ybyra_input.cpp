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

#include "options/ybyra_input.hpp"

#include "tools/misc.hpp"

#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/util/container/dataframe.hpp"
#include "genesis/util/container/dataframe/reader.hpp"
#include "genesis/util/core/fs.hpp"
#include "genesis/util/core/logging.hpp"

#include <algorithm>
#include <stdexcept>

using namespace genesis::tree;
using namespace genesis::util::container;
using namespace genesis::util::core;
using namespace genesis::util::io;

// =================================================================================================
//      Local Helpers
// =================================================================================================

/**
 * @brief Read a tab-separated table with header, and check that it has the needed columns.
 */
static Dataframe read_table_( std::string const& path, std::vector<std::string> const& columns )
{
    auto reader = DataframeReader<std::string>( '\t' ).row_names_from_first_col( false );
    auto table = reader.read( from_file( path ));
    for( auto const& column : columns ) {
        if( ! table.has_col_name( column )) {
            throw std::runtime_error(
                "Table file " + path + " does not have the expected column \"" + column + "\""
            );
        }
    }
    return table;
}

// =================================================================================================
//      Setup Functions
// =================================================================================================

void YbyraInputOptions::add_ybyra_input_opts_to_app(
    CLI::App* sub,
    bool with_calls,
    std::string const& group
) {
    with_calls_ = with_calls;

    // Calls files
    CLI::Option* calls_opt = nullptr;
    if( with_calls ) {
        calls_opt = calls_input_.add_multi_file_input_opt_to_app(
            sub, "calls", "ybyra calls", "calls", "calls", false, group,
            "These are the per-sample `calls/<sample>.calls` files of a ybyra run."
        );
    }

    // Ybyra dir
    ybyra_dir_.option = sub->add_option(
        "--ybyra-dir",
        ybyra_dir_.value,
        with_calls
        ? "Output directory of a ybyra run. All sample calls files in its `calls/` directory are "
          "used, as well as the placements in `aggregate.yplace`, if present. "
          "This is an alternative to providing the files individually."
        : "Output directory of a ybyra run. The placements are read from its `aggregate.yplace`, "
          "and the samples that could not be placed from its `fail.yplace`, if present. "
          "This is an alternative to providing the placements file."
    );
    ybyra_dir_.option->check( CLI::ExistingDirectory );
    ybyra_dir_.option->group( group );
    if( calls_opt ) {
        ybyra_dir_.option->excludes( calls_opt );
    }

    // Placements file
    placements_file_.option = sub->add_option(
        "--placements-file",
        placements_file_.value,
        with_calls
        ? "The `aggregate.yplace` file of a ybyra run, containing the placement of each sample. "
          "If provided, the placement of each sample is marked in the plot."
        : "The `aggregate.yplace` file of a ybyra run, containing the placement of each sample."
    );
    placements_file_.option->check( CLI::ExistingFile );
    placements_file_.option->group( group );
    placements_file_.option->excludes( ybyra_dir_.option );

    // Damage
    if( ! with_calls ) {
        return;
    }
    exclude_damage_.option = sub->add_flag(
        "--exclude-damage",
        exclude_damage_.value,
        "Exclude calls that are flagged as potentially affected by ancient DNA damage. "
        "This needs to match the `damage_filter` setting of the ybyra run, "
        "for the scores here to be identical to the ones computed by ybyra."
    );
    exclude_damage_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::vector<YbyraInputOptions::SampleFile> const& YbyraInputOptions::sample_files() const
{
    if( ! sample_files_.empty() ) {
        return sample_files_;
    }

    // Get the files, either from the ybyra dir, or as provided by the user.
    std::vector<std::string> paths;
    if( ybyra_dir_.is_set() ) {
        auto const calls_dir = dir_normalize_path( ybyra_dir_.value ) + "calls";
        if( ! dir_exists( calls_dir )) {
            throw CLI::ValidationError(
                "--ybyra-dir", "Directory does not contain a `calls/` directory: " + ybyra_dir_.value
            );
        }
        paths = dir_list_files( calls_dir, true, ".*\\.calls$" );
    } else if( calls_input_.provided() ) {
        paths = calls_input_.file_paths();
    } else {
        throw CLI::ValidationError(
            "Input", "Either `--calls-path` or `--ybyra-dir` has to be provided."
        );
    }
    if( paths.empty() ) {
        throw CLI::ValidationError( "Input", "No calls files found." );
    }

    // Make the list of samples, using the file names as sample names.
    for( auto const& path : paths ) {
        sample_files_.push_back({ file_basename( path, { ".calls" }), path });
    }
    std::sort(
        sample_files_.begin(), sample_files_.end(),
        []( SampleFile const& lhs, SampleFile const& rhs ){
            return lhs.name < rhs.name;
        }
    );
    for( size_t i = 1; i < sample_files_.size(); ++i ) {
        if( sample_files_[i].name == sample_files_[i-1].name ) {
            throw std::runtime_error(
                "Duplicate sample name \"" + sample_files_[i].name + "\" in calls files " +
                sample_files_[i-1].path + " and " + sample_files_[i].path
            );
        }
    }
    LOG_MSG1 << "Found " << sample_files_.size() << " sample calls file"
             << ( sample_files_.size() != 1 ? "s" : "" );
    return sample_files_;
}

std::unordered_map<std::string, YbyraInputOptions::Placement> const&
YbyraInputOptions::placements() const
{
    if( placements_read_ ) {
        return placements_;
    }
    placements_read_ = true;

    auto const path = placements_file_path_();
    if( path.empty() ) {
        return placements_;
    }

    auto const table = read_table_( path, { "individual", "optplacement", "tree_score", "flag" });
    auto const& col_ind = table[ "individual"   ].as<std::string>();
    auto const& col_opt = table[ "optplacement" ].as<std::string>();
    auto const& col_scr = table[ "tree_score"   ].as<std::string>();
    auto const& col_flg = table[ "flag"         ].as<std::string>();
    for( size_t i = 0; i < table.rows(); ++i ) {
        if( placements_.count( col_ind[i] ) > 0 ) {
            throw std::runtime_error(
                "Duplicate individual \"" + col_ind[i] + "\" in placements file " + path
            );
        }

        // Ybyra uses `...` to indicate that no flags are set.
        auto const flag = ( col_flg[i] == "..." ? "" : col_flg[i] );
        placements_[ col_ind[i] ] = Placement{ col_opt[i], std::stod( col_scr[i] ), flag };
    }
    LOG_MSG1 << "Found " << placements_.size() << " sample placement"
             << ( placements_.size() != 1 ? "s" : "" ) << " in " << path;
    return placements_;
}

std::vector<std::string> YbyraInputOptions::failed_samples() const
{
    if( ! ybyra_dir_.is_set() ) {
        return {};
    }
    auto const path = dir_normalize_path( ybyra_dir_.value ) + "fail.yplace";
    if( ! file_exists( path )) {
        return {};
    }
    auto const table = read_table_( path, { "individual" });
    auto const& column = table[ "individual" ].as<std::string>();
    return std::vector<std::string>( column.begin(), column.end() );
}

std::vector<double> YbyraInputOptions::read_sample_edge_values(
    SampleFile const& sample,
    TreeTableOptions const& tree_opts
) const {
    auto const& tree = tree_opts.get_tree();
    auto const& node_name_to_edge_index = tree_opts.node_name_to_edge_index();

    // Read the table, with the columns as written by ybyra.
    auto const table = read_table_( sample.path, { "id", "parent", "state", "damage" });
    auto const& col_id     = table[ "id"     ].as<std::string>();
    auto const& col_parent = table[ "parent" ].as<std::string>();
    auto const& col_state  = table[ "state"  ].as<std::string>();
    auto const& col_damage = table[ "damage" ].as<std::string>();

    // Add up the derived and ancestral calls per edge.
    auto edge_values = std::vector<double>( tree.edge_count(), 0.0 );
    for( size_t i = 0; i < table.rows(); ++i ) {
        if( col_damage[i] != "yes" && col_damage[i] != "no" ) {
            throw std::runtime_error(
                "Invalid damage value \"" + col_damage[i] + "\" in calls file " + sample.path
            );
        }
        if( exclude_damage_.value && col_damage[i] == "yes" ) {
            continue;
        }

        // Find the edge, and check that the parent matches, to make sure that we use the same tree.
        auto const edge_it = node_name_to_edge_index.find( col_id[i] );
        if( edge_it == node_name_to_edge_index.end() ) {
            throw std::runtime_error(
                "Node \"" + col_id[i] + "\" from calls file " + sample.path + " is not in the tree. "
                "Is the tree table the same as the one used in the ybyra run?"
            );
        }
        auto const& edge = tree.edge_at( edge_it->second );
        auto const& parent_name = edge.primary_node().data<CommonNodeData>().name;
        if( parent_name != col_parent[i] ) {
            throw std::runtime_error(
                "Node \"" + col_id[i] + "\" from calls file " + sample.path + " has parent \"" +
                col_parent[i] + "\", but in the tree, its parent is \"" + parent_name + "\". "
                "Is the tree table the same as the one used in the ybyra run?"
            );
        }

        if( col_state[i] == "derived" ) {
            edge_values[ edge_it->second ] += 1.0;
        } else if( col_state[i] == "ancestral" ) {
            edge_values[ edge_it->second ] -= 1.0;
        } else {
            throw std::runtime_error(
                "Invalid state value \"" + col_state[i] + "\" in calls file " + sample.path
            );
        }
    }
    return edge_values;
}

std::string YbyraInputOptions::mismatch_hint() const
{
    return
        "This can happen if the ybyra run used a different damage filter setting: "
        "Currently, calls flagged as damaged are " +
        std::string( exclude_damage_.value ? "excluded" : "included" ) +
        ", which can be changed with `--exclude-damage`. "
        "Alternatively, check that the tree table is the same as the one used in the ybyra run."
    ;
}

// =================================================================================================
//      Internal Helpers
// =================================================================================================

std::string YbyraInputOptions::placements_file_path_() const
{
    if( placements_file_.is_set() ) {
        return placements_file_.value;
    }
    if( ybyra_dir_.is_set() ) {
        auto const path = dir_normalize_path( ybyra_dir_.value ) + "aggregate.yplace";
        if( file_exists( path )) {
            return path;
        }
        if( ! with_calls_ ) {
            throw CLI::ValidationError(
                "--ybyra-dir", "Directory does not contain an `aggregate.yplace` file: " +
                ybyra_dir_.value
            );
        }
        LOG_WARN << "No aggregate.yplace found in " << ybyra_dir_.value
                 << ", so no placements are shown.";
    } else if( ! with_calls_ ) {
        throw CLI::ValidationError(
            "Input", "Either `--placements-file` or `--ybyra-dir` has to be provided."
        );
    }
    return "";
}
