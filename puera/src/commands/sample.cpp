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

#include "commands/sample.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/tree_drawing.hpp"

#include "genesis/tree/common_tree/tree.hpp"
// accumulate.hpp needs is_root() from function.hpp, but does not include it itself.
#include "genesis/tree/function/function.hpp"
#include "genesis/tree/function/accumulate.hpp"
#include "genesis/util/core/logging.hpp"
#include "genesis/util/text/string.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <unordered_set>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_sample( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SampleOptions>();
    auto sub = app.add_subcommand(
        "sample",
        "Plot the tree scores and placement of individual samples on the tree."
    );

    // Input
    options->tree_table.add_tree_table_opt_to_app( sub );
    options->snps_table.add_snps_table_opt_to_app( sub, false );
    options->ybyra_input.add_ybyra_input_opts_to_app( sub );

    // Drawing
    options->svg_tree.add_svg_tree_output_opts_to_app( sub );
    options->annotation.add_clade_labels_opts_to_app( sub );
    options->annotation.add_marker_opts_to_app( sub );
    options->annotation.add_title_opt_to_app( sub );

    // Color
    options->color_map.add_color_list_opt_to_app( sub, "viridis" );
    options->color_map.add_under_opt_to_app( sub, "#B40041" );
    options->color_map.add_over_opt_to_app( sub );
    options->color_norm.add_min_value_opt_to_app( sub );
    options->color_norm.add_max_value_opt_to_app( sub );
    options->color_range.option = sub->add_option(
        "--color-range",
        options->color_range.value,
        "Determine the maximum of the color scale from the scores of each sample individually "
        "(`per-sample`), or from the highest score across all samples (`all-samples`), "
        "for colors that are comparable between samples. Ignored if `--max-value` is set. "
        "Scores below the minimum, by default 0, are drawn in the `--under-color`."
    );
    options->color_range.option->group( "Color" );
    options->color_range.option->transform(
        CLI::IsMember({ "per-sample", "all-samples" }, CLI::ignore_case )
    );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( puera_cli_callback(
        sub,
        [ options ]() {
            run_sample( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

/**
 * @brief Compute the tree score of each node of the tree, indexed by the edge leading to the node.
 *
 * This is the sum of derived (+1) and ancestral (-1) calls on the path from the root to the node,
 * which is the same as the `tree_score` of ybyra.
 */
static std::vector<double> compute_tree_scores(
    SampleOptions const& options,
    YbyraInputOptions::SampleFile const& sample
) {
    auto const edge_values = options.ybyra_input.read_sample_edge_values(
        sample, options.tree_table
    );
    return genesis::tree::accumulate_edge_values_outwards(
        options.tree_table.get_tree(), edge_values
    );
}

static std::string format_score( double score )
{
    return std::to_string( static_cast<long long>( std::llround( score )));
}

/**
 * @brief Warn about samples that have a placement, but no calls file.
 */
static void check_placements_without_calls( SampleOptions const& options )
{
    std::unordered_set<std::string> sample_names;
    for( auto const& sample : options.ybyra_input.sample_files() ) {
        sample_names.insert( sample.name );
    }
    std::vector<std::string> missing;
    for( auto const& placement : options.ybyra_input.placements() ) {
        if( sample_names.count( placement.first ) == 0 ) {
            missing.push_back( placement.first );
        }
    }
    if( ! missing.empty() ) {
        std::sort( missing.begin(), missing.end() );
        LOG_WARN << "Warning: " << missing.size() << " sample"
                 << ( missing.size() != 1 ? "s have" : " has" )
                 << " a placement, but no calls file, and will not be plotted: "
                 << genesis::util::text::join( missing, ", " );
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_sample( SampleOptions const& options )
{
    using namespace genesis::util::format;

    // Get the samples and their placements, and check that we can write the output files.
    auto const& samples = options.ybyra_input.sample_files();
    auto const& placements = options.ybyra_input.placements();
    std::vector<std::string> infixes;
    for( auto const& sample : samples ) {
        infixes.push_back( "sample-" + sample.name );
    }
    options.file_output.check_output_files_nonexistence( infixes, "svg" );
    check_placements_without_calls( options );

    // Prepare everything that is the same for all samples.
    auto const setup = prepare_tree_drawing(
        options.tree_table, options.snps_table, options.svg_tree, options.annotation
    );
    auto const& tree = setup.tree;
    auto const& color_map = options.color_map.color_map();

    // If requested, get the highest score across all samples, for the maximum of the color scale.
    double all_samples_max = 0.0;
    if( options.color_range.value == "all-samples" && ! options.color_norm.max_value_is_set() ) {
        for( auto const& sample : samples ) {
            auto const scores = compute_tree_scores( options, sample );
            all_samples_max = std::max(
                all_samples_max, *std::max_element( scores.begin(), scores.end() )
            );
        }
    }

    // Plot each sample.
    for( auto const& sample : samples ) {
        LOG_MSG1 << "Plotting sample " << sample.name;
        auto const scores = compute_tree_scores( options, sample );
        auto const max_score = *std::max_element( scores.begin(), scores.end() );
        auto const min_score = *std::min_element( scores.begin(), scores.end() );
        if( max_score == 0.0 && min_score == 0.0 ) {
            LOG_WARN << "Warning: Sample " << sample.name << " has no calls, and is not plotted.";
            continue;
        }

        // Mark the placement of the sample, if there is one, and make the title.
        std::string title = sample.name;
        auto node_shapes = std::vector<SvgGroup>( tree.node_count() );
        auto const placement_it = placements.find( sample.name );
        if( placement_it != placements.end() ) {
            auto const& placement = placement_it->second;
            auto const edge_index = placement_edge_index(
                options.tree_table, placement.node, sample.name
            );

            // Check that our score at the placement node is the one that ybyra computed.
            auto const score = scores[ edge_index ];
            if( score != placement.tree_score ) {
                LOG_WARN << "Warning: The tree score of sample " << sample.name
                         << " at its placement " << placement.node << " is " << format_score( score )
                         << ", but ybyra computed " << format_score( placement.tree_score ) << ". "
                         << options.ybyra_input.mismatch_hint();
            }

            // Ties are resolved by ybyra by placing the sample at the most recent common ancestor
            // of the tied nodes, which we indicate with a dashed ring.
            bool const dashed = placement.flag.find( "most_recent_common_parent" ) != std::string::npos;
            auto const node_index = tree.edge_at( edge_index ).secondary_node().index();
            node_shapes[ node_index ] = options.annotation.make_placement_marker( setup.size_unit, dashed );

            title += " · " + placement.node + " · score " + format_score( placement.tree_score );
            if( ! placement.flag.empty() ) {
                title += " · " + placement.flag;
            }
        } else if( ! placements.empty() ) {
            title += " · not placed";
        }

        // Make the color normalization. By default, scores below zero get the under color.
        auto norm = options.color_norm.get_sequential_norm();
        if( ! options.color_norm.min_value_is_set() ) {
            norm->min_value( 0.0 );
        }
        if( ! options.color_norm.max_value_is_set() ) {
            norm->max_value(
                options.color_range.value == "all-samples" ? all_samples_max : max_score
            );
        }
        if( norm->max_value() <= norm->min_value() ) {
            norm->max_value( norm->min_value() + 1.0 );
        }

        // Draw the tree and write it.
        write_tree_svg(
            setup, scores, color_map, *norm, node_shapes,
            options.annotation, title, options.file_output, "sample-" + sample.name
        );
    }
}
