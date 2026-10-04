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

#include "commands/summary.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/tree_drawing.hpp"

#include "genesis/tree/common_tree/tree.hpp"
// accumulate.hpp needs is_root() from function.hpp, but does not include it itself.
#include "genesis/tree/function/function.hpp"
#include "genesis/tree/function/accumulate.hpp"
#include "genesis/tree/iterator/preorder.hpp"
#include "genesis/util/container/dataframe.hpp"
#include "genesis/util/container/dataframe/reader.hpp"
#include "genesis/util/core/logging.hpp"
#include "genesis/util/io/input_source.hpp"
#include "genesis/util/text/string.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_summary( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SummaryOptions>();
    auto sub = app.add_subcommand(
        "summary",
        "Plot a summary of the placements of all samples on the tree."
    );

    // Input
    options->tree_table.add_tree_table_opt_to_app( sub );
    options->snps_table.add_snps_table_opt_to_app( sub, false );
    options->ybyra_input.add_ybyra_input_opts_to_app( sub, false );
    options->groups_file.option = sub->add_option(
        "--groups-file",
        options->groups_file.value,
        "Tab-separated table with columns `sample` and `group`, with a header line. "
        "If provided, an additional summary is plotted for the samples of each group."
    );
    options->groups_file.option->check( CLI::ExistingFile );
    options->groups_file.option->group( "Input" );
    options->exclude_flags.option = sub->add_option(
        "--exclude-flags",
        options->exclude_flags.value,
        "Comma-separated list of ybyra flags, such as `low_tree_score`. "
        "Samples whose placement has any of these flags are not counted."
    );
    options->exclude_flags.option->delimiter( ',' );
    options->exclude_flags.option->group( "Input" );

    // Counting
    options->count_mode.option = sub->add_option(
        "--count-mode",
        options->count_mode.value,
        "Color each branch by the number of samples placed in the whole clade below it "
        "(`cumulative`), which shows the paths from the root to the populated clades, "
        "or by the number of samples placed at the node of the branch itself (`per-edge`)."
    );
    options->count_mode.option->group( "Counting" );
    options->count_mode.option->transform(
        CLI::IsMember({ "cumulative", "per-edge" }, CLI::ignore_case )
    );

    // Drawing
    options->svg_tree.add_svg_tree_output_opts_to_app( sub );
    options->annotation.add_clade_labels_opts_to_app( sub );
    options->count_circles.option = sub->add_flag(
        "--count-circles",
        options->count_circles.value,
        "Draw a circle at each node with placed samples, with its area proportional to the "
        "number of samples placed at that node, independently of the `--count-mode`."
    );
    options->count_circles.option->group( "Annotation" );
    options->annotation.add_marker_opts_to_app( sub );
    options->annotation.add_title_opt_to_app( sub );

    // Color
    options->color_map.add_color_list_opt_to_app( sub, "bupubk" );
    options->color_map.add_under_opt_to_app( sub, "#D0D0D0" );
    options->color_map.add_over_opt_to_app( sub );
    options->color_norm.add_scaling_opt_to_app(
        sub,
        "With `auto`, the scale is logarithmic for the `cumulative` count mode, where counts "
        "near the root are much higher than near the tips, and linear for `per-edge`."
    );
    options->color_norm.add_min_value_opt_to_app( sub );
    options->color_norm.add_max_value_opt_to_app( sub );
    options->color_range.option = sub->add_option(
        "--color-range",
        options->color_range.value,
        "With `--groups-file`, determine the maximum of the color scale and the circle sizes "
        "of each group plot from the highest count across all groups (`all-groups`), "
        "for plots that are comparable between groups, or from each group individually "
        "(`per-group`). The plot of all samples always uses its own maximum. "
        "Branches with counts below the minimum, by default 1, are drawn in the `--under-color`."
    );
    options->color_range.option->group( "Color" );
    options->color_range.option->transform(
        CLI::IsMember({ "all-groups", "per-group" }, CLI::ignore_case )
    );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->write_table.option = sub->add_flag(
        "--write-table",
        options->write_table.value,
        "Also write a table with the counts per node, for all nodes with samples placed "
        "in their clade, in tree preorder."
    );
    options->write_table.option->group( "Output" );
    options->table_format.option = sub->add_option(
        "--table-format",
        options->table_format.value,
        "Format of the `--write-table` table. With `wide`, there is one row per node, and the "
        "counts per group are additional columns. With `long`, there is one row per node and "
        "group, with group `all` for all samples."
    );
    options->table_format.option->group( "Output" );
    options->table_format.option->transform(
        CLI::IsMember({ "wide", "long" }, CLI::ignore_case )
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( puera_cli_callback(
        sub,
        [ options ]() {
            run_summary( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

/**
 * @brief Counts of the samples of one plot: all samples, or the samples of one group.
 */
struct SampleSet
{
    std::string name;
    std::vector<double> placed;
    std::vector<double> cumulative;
    size_t counted  = 0;
    size_t excluded = 0;
    size_t failed   = 0;
};

/**
 * @brief Sample groups, in the order of their first appearance in the groups file.
 */
struct Groups
{
    std::vector<std::string> names;
    std::unordered_map<std::string, std::string> sample_to_group;
};

static Groups read_groups( SummaryOptions const& options )
{
    using namespace genesis::util::container;
    using namespace genesis::util::io;

    Groups result;
    if( ! options.groups_file.is_set() ) {
        return result;
    }
    auto const& path = options.groups_file.value;

    auto reader = DataframeReader<std::string>( '\t' ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( path ));
    if( ! table.has_col_name( "sample" ) || ! table.has_col_name( "group" )) {
        throw std::runtime_error(
            "Groups file " + path + " needs a header line with columns `sample` and `group`"
        );
    }
    auto const& col_sample = table[ "sample" ].as<std::string>();
    auto const& col_group  = table[ "group"  ].as<std::string>();

    for( size_t i = 0; i < table.rows(); ++i ) {
        auto const& sample = col_sample[i];
        auto const& group  = col_group[i];
        if( result.sample_to_group.count( sample ) > 0 ) {
            throw std::runtime_error(
                "Sample \"" + sample + "\" is listed multiple times in groups file " + path
            );
        }

        // Group names are used in file names and as the group of all samples in the table.
        bool const valid = ! group.empty() && std::all_of(
            group.begin(), group.end(), []( char c ){
                return std::isalnum( static_cast<unsigned char>( c )) ||
                       c == '.' || c == '_' || c == '-';
            }
        );
        if( ! valid || group == "all" ) {
            throw std::runtime_error(
                "Invalid group name \"" + group + "\" in groups file " + path + ". "
                "Group names can only contain letters, digits, `.`, `_`, and `-`, "
                "and `all` is reserved for all samples."
            );
        }

        if( std::find( result.names.begin(), result.names.end(), group ) == result.names.end() ) {
            result.names.push_back( group );
        }
        result.sample_to_group[ sample ] = group;
    }
    LOG_MSG1 << "Found " << result.sample_to_group.size() << " samples in "
             << result.names.size() << " group" << ( result.names.size() != 1 ? "s" : "" );
    return result;
}

static std::string count_str( size_t count, std::string const& noun )
{
    return std::to_string( count ) + " " + noun + ( count != 1 ? "s" : "" );
}

/**
 * @brief Warn about samples, listing the first few of their names, sorted.
 */
static void warn_samples( std::vector<std::string> names, std::string const& message )
{
    size_t const max_names = 10;
    if( names.empty() ) {
        return;
    }
    auto const count = names.size();
    std::sort( names.begin(), names.end() );
    names.resize( std::min( count, max_names ));
    LOG_WARN << "Warning: " << count_str( count, "sample" ) << " " << message << ": "
             << genesis::util::text::join( names, ", " )
             << ( count > max_names ? ", and " + std::to_string( count - max_names ) + " more" : "" );
}

static std::string format_count( double value )
{
    return std::to_string( static_cast<long long>( std::llround( value )));
}

static void write_table( SummaryOptions const& options, std::vector<SampleSet> const& sets )
{
    using namespace genesis::tree;

    auto const& tree = options.tree_table.get_tree();
    auto target = options.file_output.get_output_target( "summary", "tsv" );
    auto& os = target->ostream();
    bool const wide = options.table_format.value == "wide";

    // Header
    if( wide ) {
        os << "node\tplaced\tcumulative";
        for( size_t i = 1; i < sets.size(); ++i ) {
            os << "\tplaced:" << sets[i].name << "\tcumulative:" << sets[i].name;
        }
        os << "\n";
    } else {
        os << "node\tgroup\tplaced\tcumulative\n";
    }

    // Rows, for all nodes that have samples placed in their clade, which are the same
    // for all groups, as the groups are subsets of all samples.
    for( auto const& it : preorder( tree )) {
        if( it.is_first_iteration() ) {
            continue;
        }
        auto const e = it.edge().index();
        if( sets[0].cumulative[e] == 0.0 ) {
            continue;
        }
        auto const& node_name = it.node().data<CommonNodeData>().name;
        if( wide ) {
            os << node_name;
            for( auto const& set : sets ) {
                os << "\t" << format_count( set.placed[e] ) << "\t" << format_count( set.cumulative[e] );
            }
            os << "\n";
        } else {
            for( auto const& set : sets ) {
                if( set.cumulative[e] == 0.0 ) {
                    continue;
                }
                os << node_name << "\t" << ( set.name.empty() ? "all" : set.name ) << "\t"
                   << format_count( set.placed[e] ) << "\t" << format_count( set.cumulative[e] )
                   << "\n";
            }
        }
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_summary( SummaryOptions const& options )
{
    using namespace genesis::tree;
    using namespace genesis::util::format;

    // Get the input, and check that we can write the output files.
    auto const& placements = options.ybyra_input.placements();
    auto const failed = options.ybyra_input.failed_samples();
    auto const groups = read_groups( options );
    std::vector<std::string> infixes = { "summary" };
    for( auto const& group : groups.names ) {
        infixes.push_back( "summary-" + group );
    }
    options.file_output.check_output_files_nonexistence( infixes, "svg" );
    if( options.write_table.value ) {
        options.file_output.check_output_files_nonexistence( "summary", "tsv" );
    }

    // Prepare the tree drawing.
    auto const setup = prepare_tree_drawing(
        options.tree_table, options.snps_table, options.svg_tree, options.annotation
    );
    auto const& tree = setup.tree;
    auto const& root_name = tree.root_node().data<CommonNodeData>().name;

    // Make the sets of samples to plot: All samples first, then each group.
    std::vector<SampleSet> sets( 1 + groups.names.size() );
    std::unordered_map<std::string, size_t> group_to_set_index;
    for( size_t i = 0; i < sets.size(); ++i ) {
        if( i > 0 ) {
            sets[i].name = groups.names[i-1];
            group_to_set_index[ sets[i].name ] = i;
        }
        sets[i].placed = std::vector<double>( tree.edge_count(), 0.0 );
    }

    // Get the indices of the sets that a sample belongs to.
    std::vector<std::string> ungrouped;
    auto get_set_indices = [&]( std::string const& sample ){
        std::vector<size_t> result = { 0 };
        if( groups.names.empty() ) {
            return result;
        }
        auto const it = groups.sample_to_group.find( sample );
        if( it == groups.sample_to_group.end() ) {
            ungrouped.push_back( sample );
        } else {
            result.push_back( group_to_set_index.at( it->second ));
        }
        return result;
    };

    // Count the placed samples per edge, except for the excluded ones.
    std::unordered_set<std::string> const excluded_flags(
        options.exclude_flags.value.begin(), options.exclude_flags.value.end()
    );
    std::unordered_set<std::string> seen_flags;
    std::vector<std::string> root_placed;
    for( auto const& entry : placements ) {
        auto const& sample = entry.first;
        auto const& placement = entry.second;

        bool excluded = false;
        for( auto const& flag : genesis::util::text::split( placement.flag, ";", true )) {
            seen_flags.insert( flag );
            excluded |= excluded_flags.count( flag ) > 0;
        }

        // The root has no branch to color, so we cannot show samples placed there.
        size_t edge_index = 0;
        if( placement.node == root_name ) {
            root_placed.push_back( sample );
            excluded = true;
        } else {
            edge_index = placement_edge_index( options.tree_table, placement.node, sample );
        }

        for( auto const i : get_set_indices( sample )) {
            if( excluded ) {
                ++sets[i].excluded;
            } else {
                sets[i].placed[ edge_index ] += 1.0;
                ++sets[i].counted;
            }
        }
    }
    for( auto const& sample : failed ) {
        for( auto const i : get_set_indices( sample )) {
            ++sets[i].failed;
        }
    }

    // Warn about everything that might be a mistake.
    for( auto const& flag : options.exclude_flags.value ) {
        if( seen_flags.count( flag ) == 0 ) {
            LOG_WARN << "Warning: Flag \"" << flag << "\" from `--exclude-flags` does not occur "
                     << "in any of the placements. Is it spelled correctly?";
        }
    }
    warn_samples( root_placed, "placed at the root of the tree, and not counted" );
    warn_samples( ungrouped, "not listed in the groups file, and only counted for all samples" );
    std::vector<std::string> unknown;
    for( auto const& entry : groups.sample_to_group ) {
        if(
            placements.count( entry.first ) == 0 &&
            std::find( failed.begin(), failed.end(), entry.first ) == failed.end()
        ) {
            unknown.push_back( entry.first );
        }
    }
    warn_samples( unknown, "from the groups file without a placement" );

    // Accumulate the counts of each clade.
    for( auto& set : sets ) {
        set.cumulative = accumulate_edge_values_inwards( tree, set.placed );
    }

    // Get the maxima across groups, for comparable group plots.
    bool const cumulative = options.count_mode.value == "cumulative";
    auto max_of = []( std::vector<double> const& values ){
        return *std::max_element( values.begin(), values.end() );
    };
    double groups_max_value  = 0.0;
    double groups_max_placed = 0.0;
    for( size_t i = 1; i < sets.size(); ++i ) {
        auto const& values = cumulative ? sets[i].cumulative : sets[i].placed;
        groups_max_value  = std::max( groups_max_value,  max_of( values ));
        groups_max_placed = std::max( groups_max_placed, max_of( sets[i].placed ));
    }

    // Plot each set.
    auto const& color_map = options.color_map.color_map();
    for( size_t i = 0; i < sets.size(); ++i ) {
        auto const& set = sets[i];
        auto const infix = infixes[i];
        if( set.counted == 0 ) {
            LOG_WARN << "Warning: No samples counted for " << infix << ", which is not plotted.";
            continue;
        }
        LOG_MSG1 << "Plotting " << infix << " with " << count_str( set.counted, "sample" );

        auto const& values = cumulative ? set.cumulative : set.placed;
        bool const use_groups_max = i > 0 && options.color_range.value == "all-groups";
        auto const max_value  = use_groups_max ? groups_max_value  : max_of( values );
        auto const max_placed = use_groups_max ? groups_max_placed : max_of( set.placed );

        // Make the color normalization. By default, branches without samples get the under color.
        auto norm = options.color_norm.get_sequential_norm( cumulative );
        if( ! options.color_norm.min_value_is_set() ) {
            norm->min_value( 1.0 );
        }
        if( ! options.color_norm.max_value_is_set() ) {
            norm->max_value( max_value );
        }
        if( norm->max_value() <= norm->min_value() ) {
            norm->max_value( norm->min_value() + 1.0 );
        }

        // Make the circles at the nodes with placed samples.
        auto node_shapes = std::vector<SvgGroup>( tree.node_count() );
        if( options.count_circles.value ) {
            for( size_t e = 0; e < tree.edge_count(); ++e ) {
                if( set.placed[e] > 0.0 ) {
                    auto const node_index = tree.edge_at( e ).secondary_node().index();
                    node_shapes[ node_index ] = options.annotation.make_count_circle(
                        setup.size_unit, set.placed[e] / max_placed
                    );
                }
            }
        }

        // Make the title.
        std::string title = "summary";
        if( ! set.name.empty() ) {
            title += " · " + set.name;
        }
        title += " · " + count_str( set.counted, "sample" );
        if( set.excluded > 0 ) {
            title += " · " + std::to_string( set.excluded ) + " excluded";
        }
        if( set.failed > 0 ) {
            title += " · " + std::to_string( set.failed ) + " failed";
        }
        if( options.count_circles.value ) {
            title += " · circle area ∝ samples placed, max " + format_count( max_placed );
        }

        write_tree_svg(
            setup, values, color_map, *norm, node_shapes,
            options.annotation, title, options.file_output, infix
        );
    }

    if( options.write_table.value ) {
        write_table( options, sets );
    }
}
