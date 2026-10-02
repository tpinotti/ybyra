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

#include "functions/functions.hpp"

#include "genesis/tree/common_tree/tree.hpp"
#include "genesis/tree/function/functions.hpp"

using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//      Helper Functions
// =================================================================================================

std::vector<double> propagate_edge_values_inwards(
    Tree const& tree,
    std::vector<double> const& edge_values
) {
    // LOG_DBG << "propagate_edge_values_inwards";
    return accumulate_edge_values_inwards( tree, edge_values );
    // auto result = std::vector<double>( edge_values.size(), 0.0 );
    // assert( result.size() == tree.edge_count() );
    //
    // for( auto const& node_it : postorder( tree )) {
    //     if( is_root( node_it.node() )) {
    //         continue;
    //     }
    //     result[node_it.edge().index()] += edge_values[node_it.edge().index()];
    //     for( auto const& link_it : node_links( node_it.node() )) {
    //         if( link_it.is_first_iteration() ) {
    //             continue;
    //         }
    //         result[node_it.edge().index()] += result[link_it.edge().index()];
    //     }
    // }
    // return result;
}

std::vector<double> propagate_edge_snp_counts(
    Tree const& tree,
    std::vector<double> const& edge_snp_counts
) {
    // LOG_DBG << "propagate_edge_snp_counts";
    return accumulate_edge_values_outwards( tree, edge_snp_counts );

    // auto result = std::vector<double>( edge_snp_counts.size(), 0.0 );
    // assert( result.size() == tree.edge_count() );
    //
    // for( auto const& node_it : levelorder( tree )) {
    //     if( is_root( node_it.node() )) {
    //         continue;
    //     }
    //     // LOG_DBG << "at " << node_it.edge().index() << " " << node_it.node().data<CommonNodeData>().name << " with " << result[node_it.edge().index()];
    //     // if( result[node_it.edge().index()] != 0.0 ) {
    //     // }
    //     // assert( result[node_it.edge().index()] == 0.0 );
    //     assert(
    //         node_it.edge().index() != node_it.edge().primary_node().link().edge().index() ||
    //         is_root(node_it.edge().primary_node())
    //     );
    //     assert(
    //         result[node_it.edge().index()] == result[node_it.edge().primary_node().link().edge().index()] ||
    //         is_root(node_it.edge().primary_node())
    //     );
    //     result[node_it.edge().index()] += edge_snp_counts[node_it.edge().index()];
    //     for( auto const& link_it : node_links( node_it.node() )) {
    //         if( link_it.is_first_iteration() ) {
    //             // assert( result[node_it.edge().index()] == edge_snp_counts[link_it.edge().index()] );
    //             continue;
    //         }
    //         assert( result[link_it.edge().index()] == 0.0 );
    //         result[link_it.edge().index()] = result[node_it.edge().index()];
    //         // result[node_it.edge().index()] += result[link_it.edge().index()];
    //     }
    // }
    // return result;
}

/*
void make_sample_jplace(
    Tree const& tree,
    std::vector<double> const& edge_masses,
    std::string const& sample_name,
    std::string const& out_file
) {
    // LOG_MSG << "make_sample_jplace";

    // We create a fake placement sample. The edge values we want to use here are the accumulated
    // SNP counts, which we could normalize here first to obtain a single placement with LWRs.
    // However, in the final output, we probably rather want to see the SNP counts directly,
    // so instead we misuse the multiplicity for the SNP count, and just place each location
    // with an LWR of 1.
    auto sample = Sample( convert_common_tree_to_placement_tree( tree ));
    for( size_t i = 0; i < edge_masses.size(); ++i ) {
        if( !std::isfinite( edge_masses[i] ) || edge_masses[i] <= 0.0 ) {
            continue;
        }

        auto& edge = sample.tree().edge_at( i );
        auto pquery = Pquery();
        auto& placement = pquery.add_placement( edge );
        // auto placement = PqueryPlacement( edge );
        placement.like_weight_ratio = 1.0;
        placement.likelihood = 0.0;
        placement.proximal_length = edge.data<PlacementEdgeData>().branch_length / 2.0;
        placement.pendant_length = 0.0;
        auto name = PqueryName(
            sample_name + "-" + std::to_string( edge.data<PlacementEdgeData>().edge_num() ),
            edge_masses[i]
        );
        pquery.add_name( name );
        sample.add( pquery );
    }
    JplaceWriter().write( sample, to_file( out_file ));
}
*/
