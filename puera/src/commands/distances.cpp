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

#include "commands/distances.hpp"
#include "functions/functions.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

// =================================================================================================
//      Enum Mapping
// =================================================================================================

// std::vector<std::pair<std::string, genesis::utils::HeatmapParameters::ColorNorm>> const color_norm_map = {
//     { "linear",      genesis::utils::HeatmapParameters::ColorNorm::kLinear },
//     { "logarithmic", genesis::utils::HeatmapParameters::ColorNorm::kLogarithmic },
//     // { "diverging",   genesis::utils::HeatmapParameters::ColorNorm::kDiverging },
// };

// =================================================================================================
//      Setup
// =================================================================================================

void setup_distances( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<DistancesOptions>();
    auto sub = app.add_subcommand(
        "distances",
        "Create a tree plot visualizing the placement of samples."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    options->tree_table_input.add_tree_table_opt_to_app( sub );
    options->snps_table_input.add_snps_table_opt_to_app( sub );

    // // Add the file input options, separate for both types
    // options->json_input.add_multi_file_input_opt_to_app( sub, "json", "json", "json", "json", false );
    // options->csv_input.add_multi_file_input_opt_to_app( sub, "csv", "csv", "csv", "csv", false );
    // options->json_input.option()->excludes( options->csv_input.option() );
    // options->csv_input.option()->excludes( options->json_input.option() );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Output. Here, compression does not really make sense, as we want to be able to see the
    // picture files in the end anyway, so we do not add compression here.
    options->file_output.add_default_output_opts_to_app( sub );
    // options->file_output.add_file_compress_opt_to_app( sub );

    // -------------------------------------------------------------------------
    //     Callback
    // -------------------------------------------------------------------------

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( puera_cli_callback(
        sub,
        [ options ]() {
            run_distances( *options );
        }
    ));
}

// =================================================================================================
//     Tree Helper Functions
// =================================================================================================

template<typename ElementType>
void compute_leaf_snp_sums(
    TreeNode const& start_node,
    ElementType const& element,
    // std::vector<size_t> const& edge_snp_counts,
    utils::Matrix<double> const& pairwise_dists,
    std::ofstream& os
) {
    size_t leaf_cnt = 0;
    double snp_sum = 0;
    for( auto const& elem : preorder( element )) {
        if( ! is_leaf( elem.node() )) {
            continue;
        }
        ++leaf_cnt;
        snp_sum += pairwise_dists( start_node.index(), elem.node().index() );
    }
    double const avg = snp_sum / static_cast<double>( leaf_cnt );
    os << avg << "\t" << leaf_cnt << "\n";
}

void make_node_distances_tables(
    Tree const& tree,
    // std::vector<size_t> const& edge_snp_counts,
    std::string const& out_dir
) {
    // Mutation rate from https://doi.org/10.1038/ng.3171
    // 8.71 × 10−10 mutations per position per year.
    // double const mu = 8.71e-10;

    // Distance from root for every node in the tree.
    LOG_INFO << "make_node_distances_tables";
    LOG_INFO << "dists";
    auto const dists_from_root = node_branch_length_distance_vector( tree );
    auto const pairwise_dists = node_branch_length_distance_matrix( tree );

    LOG_INFO << "table";
    std::ofstream os;
    file_output_stream( out_dir + "node_snp_distances.tsv", os );
    os << "node\tdist_to_root\tavg_dist_to_leaves\tnum_leaves\n";
    for( auto const& node : tree.nodes() ) {
        os << node.data<CommonNodeData>().name << "\t";
        os << dists_from_root[ node.index() ] << "\t";
        if( is_root( node )) {
            compute_leaf_snp_sums( node, tree, pairwise_dists, os );
        } else {
            compute_leaf_snp_sums( node, Subtree( node ), pairwise_dists, os );
        }
    }
    LOG_INFO << "done";
}

// =================================================================================================
//      Run
// =================================================================================================

void run_distances( DistancesOptions const& options )
{
    using namespace genesis::tree;
    using namespace genesis::utils;

    // Check the out file target pattern, before we open any files.
    options.file_output.check_output_files_nonexistence(
        "tree", std::vector<std::string>{ "bmp", "svg" }
    );


    // Get the tree
    auto const& tree = options.tree_table_input.get_tree();
    // read_tree_from_table_file( tree_table_infile );

    // CommonTreeNewickWriter().write( tree, utils::to_file( tree_table_infile + ".newick" ));


    // Set the branch lengths to the numbers of snps on each of them
    auto const edge_snp_counts = options.snps_table_input.get_snps_per_edge_counts(
        options.tree_table_input
    );
    options.tree_table_input.set_tree_branch_length_to_snp_counts( edge_snp_counts );

   // Make a pairwise distance table for the tree.
    make_node_distances_tables( tree, out_dir );
}
