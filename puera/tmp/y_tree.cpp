/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2024 Lucas Czech

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

#include "genesis/genesis.hpp"

#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     Tree Table
// =================================================================================================

// Custom hash function for std::pair<std::string, std::string>
struct pair_hash {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const {
        auto hash1 = std::hash<std::string>{}(p.first);
        auto hash2 = std::hash<std::string>{}(p.second);
        return hash1 ^ (hash2 << 1); // Combine the two hash values
    }
};

// Custom equality function for std::pair<std::string, std::string>
struct pair_equal {
    bool operator()(const std::pair<std::string, std::string>& p1,
                    const std::pair<std::string, std::string>& p2) const {
        return p1.first == p2.first && p1.second == p2.second;
    }
};

Tree read_tree_from_table_file(
    std::string const& table_file
) {
    LOG_DBG << "read_tree_from_table_file";

    // Read the input table into a dataframe.
    auto reader = DataframeReader<std::string>( ',' ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( table_file ));
    LOG_DBG << "tree table columns: " << join( table.col_names() );
    auto const& raw_children = table["id"].as<std::string>().to_vector();
    auto const& raw_parents = table["parent"].as<std::string>().to_vector();

    return make_tree_from_parents_table( raw_children, raw_parents );

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

std::unordered_map<std::string, size_t> make_node_name_to_edge_index( Tree const& tree )
{
    // Create a map of branch names to edge indices, for speed
    std::unordered_map<std::string, size_t> result;
    for( auto const& node : tree.nodes() ) {
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
//     SNP Table
// =================================================================================================

std::vector<size_t> get_snp_to_edge_counts(
    Tree const& tree,
    std::string const& table_file
) {
    LOG_DBG << "get_snp_to_edge_counts";

    // Create a map of branch names to edge indices, for speed
    auto const node_name_to_edge_index = make_node_name_to_edge_index( tree );

    // First read the table containing the snps and branches
    auto reader = DataframeReader<std::string>( '\t' ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( table_file ));
    LOG_DBG << "snp table columns: " << join( table.col_names() );

    // Shortcuts for the columsn of the table that we need.
    auto const& col_pos    = table[ "position" ].as<std::string>().to_vector();
    auto const& col_id     = table[ "id"       ].as<std::string>().to_vector();
    auto const& col_parent = table[ "parent"   ].as<std::string>().to_vector();
    if( col_pos.size() != col_id.size() || col_pos.size() != col_parent.size() ) {
        throw std::runtime_error( "Wrong column sizes" );
    }

    // Finally, create a vector of edge indices counting the snps on each edge.
    auto edge_snp_counts = std::vector<size_t>( tree.edge_count() );
    size_t cnt = 0;
    for( size_t i = 0; i < col_pos.size(); ++i ) {
        if( col_id[i].empty() ) {
            LOG_WARN << "empty snp table id at " << i;
            continue;
        }

        // LOG_DBG << col_pos[i];
        if( node_name_to_edge_index.count( col_id[i] ) == 0 ) {
            throw std::runtime_error( "No child with name " + col_id[i] + " at snpID " + col_pos[i] );
        }

        // auto const pos = std::stoul( col_pos[i] );
        auto const edge_index = node_name_to_edge_index.at( col_id[i] );

        auto const& parent_name = tree.edge_at( edge_index ).primary_node().data<CommonNodeData>().name;
        auto const& child_name = tree.edge_at( edge_index ).secondary_node().data<CommonNodeData>().name;
        if( parent_name != col_parent[i] ) {
            throw std::runtime_error(
                "Wrong parent name " + parent_name + " instead of " + col_parent[i] +
                " at snpID " + col_pos[i]
            );
        }
        if( child_name != col_id[i] ) {
            throw std::runtime_error(
                "Wrong child name " + child_name + " instead of " + col_id[i] +
                " at snpID " + col_pos[i]
            );
        }
        ++edge_snp_counts[ edge_index ];
        ++cnt;
    }
    LOG_INFO << "used " << cnt << " entries in snp table";
    return edge_snp_counts;
}

// =================================================================================================
//     SNPs on Tree Helper Functions
// =================================================================================================

std::unordered_map<size_t, size_t> read_snp_to_edge_map(
    Tree const& tree,
    std::string const& table_file
) {
    LOG_DBG << "read_snp_to_edge_map";

    // Create a map of branch names to edge indices, for speed
    auto const node_name_to_edge_index = make_node_name_to_edge_index( tree );

    // First read the table containing the snps and branches
    auto reader = DataframeReader<std::string>( '\t' ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( table_file ));
    // LOG_DBG << "snp table columns: " << join( table.col_names() );

    // Shortcuts for the columsn of the table that we need.
    auto const& col_pos    = table[ "position" ].as<std::string>().to_vector();
    auto const& col_id     = table[ "id"       ].as<std::string>().to_vector();
    auto const& col_parent = table[ "parent"   ].as<std::string>().to_vector();
    if( col_pos.size() != col_id.size() || col_pos.size() != col_parent.size() ) {
        throw std::runtime_error( "Wrong column sizes" );
    }

    // Finally, create the map from snp positions to the branch they appear on.
    std::unordered_map<size_t, size_t> snp_to_edge_index;
    for( size_t i = 0; i < col_pos.size(); ++i ) {
        // LOG_DBG << col_pos[i];
        if( node_name_to_edge_index.count( col_id[i] ) == 0 ) {
            throw std::runtime_error( "No child with name " + col_id[i] + " at snpID " + col_pos[i] );
        }

        auto const pos = std::stoul( col_pos[i] );
        auto const edge_index = node_name_to_edge_index.at( col_id[i] );
        if( snp_to_edge_index.count( pos ) > 0 ) {
            // throw std::runtime_error( "Duplicate position: " + col_pos[i] );
            if( snp_to_edge_index[ pos ] != edge_index ) {
                throw std::runtime_error(
                    "Duplicate position " + col_pos[i] + " with different edge indices"
                );
            }
        }

        auto const& parent_name = tree.edge_at( edge_index ).primary_node().data<CommonNodeData>().name;
        auto const& child_name = tree.edge_at( edge_index ).secondary_node().data<CommonNodeData>().name;
        if( parent_name != col_parent[i] ) {
            throw std::runtime_error(
                "Wrong parent name " + parent_name + " instead of " + col_parent[i] +
                " at snpID " + col_pos[i]
            );
        }
        if( child_name != col_id[i] ) {
            throw std::runtime_error(
                "Wrong child name " + child_name + " instead of " + col_id[i] +
                " at snpID " + col_pos[i]
            );
        }
        snp_to_edge_index[ pos ] = edge_index;
    }
    return snp_to_edge_index;
}

void set_tree_branch_length_to_snp_counts(
    std::vector<size_t> const& edge_snp_counts,
    Tree& tree
) {
    assert( edge_snp_counts.size() == tree.edge_count() );
    for( size_t i = 0; i < edge_snp_counts.size(); ++i ) {
        tree.edge_at( i ).data<CommonEdgeData>().branch_length = edge_snp_counts[i];
    }
}

std::vector<double> make_sample_edge_snp_counts(
    Tree const& tree,
    std::unordered_map<size_t, size_t> const& snp_to_edge_index,
    std::string const& table_file
) {
    LOG_DBG << "make_sample_edge_snp_counts";

    // Read the table
    auto reader = DataframeReader<std::string>( '\t' );
    reader.row_names_from_first_col( false );
    reader.col_names_from_first_row( false );
    auto const table = reader.read( from_file( table_file ));
    auto const& col_pos = table[ 0 ].as<std::string>().to_vector();

    // Prepare the result vector and loop the table
    auto edge_values = std::vector<double>( tree.edge_count() );
    for( size_t i = 0; i < table.rows(); ++i ) {
        // LOG_DBG << col_pos[i];
        auto const pos = std::stoul( col_pos[i] );
        if( snp_to_edge_index.count( pos ) == 0 ) {
            throw std::runtime_error( "Invalid position " + col_pos[i] );
        }
        edge_values[snp_to_edge_index.at(pos)] += 1.0;
    }
    return edge_values;
}

std::vector<double> make_sample_edge_snp_counts_derived(
    Tree const& tree,
    std::string const& sample_table_infile,
    bool exclude_damage
) {
    LOG_DBG << "make_sample_edge_snp_counts_derived";

    // Create a map of branch names to edge indices, for speed
    auto const node_name_to_edge_index = make_node_name_to_edge_index( tree );

    // Read the table
    auto reader = DataframeReader<std::string>( '\t' );
    reader.row_names_from_first_col( false );
    // reader.col_names_from_first_row( false );
    auto const table = reader.read( from_file( sample_table_infile ));

    // Shortcuts for the columsn of the table that we need.
    auto const& col_snpid  = table[ "snpId"  ].as<std::string>().to_vector();
    auto const& col_id     = table[ "id"     ].as<std::string>().to_vector();
    auto const& col_parent = table[ "parent" ].as<std::string>().to_vector();
    auto const& col_state  = table[ "state"  ].as<std::string>().to_vector();
    auto const& col_damage = table[ "damage" ].as<std::string>().to_vector();

    // Prepare the result vector and loop the table
    auto edge_values = std::vector<double>( tree.edge_count() );
    size_t cnt = 0;
    size_t excl = 0;
    for( size_t i = 0; i < table.rows(); ++i ) {
        if( col_id[i].empty() ) {
            LOG_WARN << "empty col_id[i] at " << i;
            continue;
        }
        if( exclude_damage && col_damage[i] == "yes" ) {
            ++excl;
            continue;
        }
        if( col_damage[i] != "yes" && col_damage[i] != "no" ) {
            LOG_WARN << "col damage: \"" << col_damage[i] << "\" at " << i;
        }
        if( node_name_to_edge_index.count( col_id[i] ) == 0 ) {
            throw std::runtime_error( "No child with name " + col_id[i] + " at snpID " + col_snpid[i] );
        }
        auto const edge_index = node_name_to_edge_index.at( col_id[i] );

        auto const& parent_name = tree.edge_at( edge_index ).primary_node().data<CommonNodeData>().name;
        auto const& child_name = tree.edge_at( edge_index ).secondary_node().data<CommonNodeData>().name;
        if( parent_name != col_parent[i] ) {
            throw std::runtime_error(
                "Wrong parent name " + parent_name + " instead of " + col_parent[i] +
                " at snpID " + col_snpid[i]
            );
        }
        if( child_name != col_id[i] ) {
            throw std::runtime_error(
                "Wrong child name " + child_name + " instead of " + col_id[i] +
                " at snpID " + col_snpid[i]
            );
        }

        if( col_state[i] == "derived" || col_state[i] == "+" ) {
            edge_values[edge_index] += 1.0;
        } else if( col_state[i] == "ancestral" || col_state[i] == "-" ) {
            edge_values[edge_index] -= 1.0;
        } else {
            LOG_WARN << "col state: \"" << col_state[i] << "\" at " << i;
        }
        ++cnt;
    }
    LOG_INFO << "used " << cnt << " rows of sample table";
    if( excl > 0 ) {
        LOG_INFO << "excluded " << excl << " rows with damage";
    }
    return edge_values;
}

std::vector<double> propagate_edge_snp_counts(
    Tree const& tree,
    std::vector<double> const& edge_snp_counts
) {
    LOG_DBG << "propagate_edge_snp_counts";
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

void write_edge_snp_values(
    Tree const& tree,
    std::vector<double> const& edge_values,
    std::vector<double> const& clade_values,
    std::string const& outfile
) {
    assert( edge_values.size() == tree.edge_count() );
    assert( clade_values.size() == tree.edge_count() );
    LOG_DBG << "write_edge_snp_values";

    std::ofstream ofs;
    utils::file_output_stream( outfile, ofs );
    ofs << "id\tedge\tsum\n";
    for( size_t i = 0; i < tree.edge_count(); ++i ) {
        auto const& child_name = tree.edge_at( i ).secondary_node().data<CommonNodeData>().name;
        ofs << child_name << "\t" << edge_values[i] << "\t" << clade_values[i] << "\n";
    }
}

// =================================================================================================
//     Tree Visualization
// =================================================================================================

std::vector<utils::SvgGroup> get_max_value_edge_marker(
    Tree const& tree,
    std::vector<double> const& clade_values,
    bool place_marker = true
) {
    assert( tree.edge_count() == clade_values.size() );

    // Get the edge index with the highest value.
    auto max_index = std::distance(
        clade_values.begin(), std::max_element( clade_values.begin(), clade_values.end() )
    );
    // size_t max_i;
    // double max_v = 0.0;
    // for( size_t i = 0; i < clade_values.size(); ++i ) {
    //     if( clade_values[i] > max_v ) {
    //         max_v = clade_values[i];
    //         max_i = i;
    //     }
    // }

    // Now max_i is the edge with the highest value.
    auto const orange = color_from_hex( "#FF6600" );
    auto edge_shapes = std::vector<utils::SvgGroup>( tree.edge_count() );
    if( place_marker ) {
        edge_shapes[ max_index ].add(
            SvgCircle( SvgPoint( 0,0 ), 500, SvgStroke( orange, 150 ), SvgFill( SvgFill::Type::kNone ))
        );
    }
    return edge_shapes;
}

std::vector<utils::SvgGroup> get_node_labels_and_remove_names(
    Tree& tree
) {
    // Now max_i is the edge with the highest value.
    auto node_shapes = std::vector<utils::SvgGroup>( tree.node_count() );
    for( auto& node : tree.nodes() ) {
        auto& name = node.data<CommonNodeData>().name;
        if( name.size() == 1 ) {
            auto text = SvgText( name );
            text.font.size = 1000;
            text.anchor = SvgText::Anchor::kMiddle;
            text.dominant_baseline = SvgText::DominantBaseline::kMiddle;
            node_shapes[node.index()].add( text );
        }
        name = "";
    }
    return node_shapes;
}

void add_edge_name_labels_and_remove_names(
    Tree& tree,
    std::vector<utils::SvgGroup>& edge_shapes
) {
    for( auto& node : tree.nodes() ) {
        auto& name = node.data<CommonNodeData>().name;
        if( name.size() == 1 ) {
            // Make a big label with the node name
            auto text = SvgText( name );
            text.font.size = 1500;
            text.anchor = SvgText::Anchor::kMiddle;
            text.dominant_baseline = SvgText::DominantBaseline::kCentral;
            text.dy = ".08em";

            // Add a white half transparent circle and the label to the edge
            auto& target_edge = edge_shapes[node.primary_link().edge().index()];
            target_edge.add( SvgCircle(
                SvgPoint( 0,0 ), 750,
                SvgStroke( SvgStroke::Type::kNone ),
                SvgFill( Color( 1.0, 1.0, 1.0, 0.7 ) )
            ));
            target_edge.add( text );
        }
        name = "";
    }
}

void mark_furthest_leaf(
    Tree& tree,
    std::vector<utils::SvgGroup>& node_shapes
) {
    // Get the node index of the node that is furthest from the root.
    auto const dists = node_branch_length_distance_vector( tree );
    auto max_index = std::distance(
        dists.begin(), std::max_element( dists.begin(), dists.end() )
    );

    // Add a red dot there, and print out the name.
    (void) node_shapes;
    // node_shapes[ max_index ].add( SvgCircle(
    //     SvgPoint( 0,0 ), 100,
    //     SvgStroke( SvgStroke::Type::kNone ),
    //     SvgFill( color_from_hex( "#B40041" ) )
    // ));
    LOG_INFO << "Furthest leaf: " << tree.node_at(max_index).data<CommonNodeData>().name;
}

struct CountSvgSettings
{
    bool with_negatives = false;
    bool phylogram = true;
    bool with_max_marker = true;
};

void sample_counts_tree_to_svg(
    Tree tree,
    std::vector<double> const& edge_values,
    std::string const& out_base,
    CountSvgSettings settings = CountSvgSettings()
) {
    LOG_DBG << "write_color_tree_to_svg_file";

    // General layout params
    LayoutParameters params;
    params.shape = LayoutShape::kCircular;
    if( settings.phylogram ) {
        params.type = LayoutType::kPhylogram;
    } else {
        params.type = LayoutType::kCladogram;
    }
    params.stroke.width = 150;

    // Get an edge marker for the highest value edge, and use nice node names
    // auto const node_markers = get_node_labels_and_remove_names( tree );
    // auto const edge_markers = get_max_value_edge_marker( tree, edge_values );
    auto node_markers = std::vector<SvgGroup>( tree.node_count() );
    auto edge_markers = get_max_value_edge_marker( tree, edge_values, settings.with_max_marker );

    // Some more manipulations on the visualization
    // mark_furthest_leaf( tree, node_markers );

    // Lastly, we remove the remaining node names, so that they don't show up in the tree,
    // only kepping nice big labels for the main branches.
    add_edge_name_labels_and_remove_names( tree, edge_markers );

    if( settings.with_negatives ) {
        auto color_map = ColorMap( color_list_spectral() );
        auto color_norm = ColorNormalizationDiverging( edge_values );
        auto color_vals = color_map( color_norm, edge_values );
        write_color_tree_to_svg_file(
            tree, params, color_vals, color_map, color_norm,
            node_markers, edge_markers, out_base + ".svg"
        );
    } else {
        // auto color_map = ColorMap( color_list_bupubk() );
        // auto color_norm = ColorNormalizationLinear( edge_values );

        // Make a color map according to the derived and ancestral counts
        auto color_map = ColorMap( color_list_viridis() );
        auto color_norm = ColorNormalizationLinear();
        color_norm.autoscale_max( edge_values );
        color_map.under_color( color_from_hex( "#B40041" ));
        auto color_vals = color_map( color_norm, edge_values );

        // Write the tree
        write_color_tree_to_svg_file(
            tree, params, color_vals, color_map, color_norm,
            node_markers, edge_markers, out_base + ".svg"
        );
    }
}

// =================================================================================================
//     Run Single Sample Visualization
// =================================================================================================

void run_single_sample(
    std::string const& tree_table_infile,
    std::string const& snps_table_infile,
    std::string const& sample_table_infile,
    bool exclude_damage
) {
    LOG_DBG << "============ run " << ( exclude_damage ? "no damage" : "all" );

    // Get the output base name
    auto const out_dir = file_path( real_path( sample_table_infile )) + "/out/";
    std::string const dmg_suffix = exclude_damage ? ".no-damage" : ".all";
    auto const out_base = out_dir + file_basename( sample_table_infile ) + dmg_suffix;
    LOG_DBG << "writing files to " << out_base;

    // Get the tree
    auto tree = read_tree_from_table_file( tree_table_infile );
    // CommonTreeNewickWriter().write( tree, utils::to_file( tree_table_infile + ".newick" ));

    // Get the mapping of snp positions to edge indices, and the counts for the sample.
    // Forgot if that was just a test or needed... will have to check later... :-(
    // auto const snp_to_edge_index = read_snp_to_edge_map( tree, snps_table_infile );
    // auto const edge_values = make_sample_edge_snp_counts(
    //     tree, snp_to_edge_index, sample_table_infile
    // );

    // Set the branch lengths to the numbers of snps on each of them
    auto const edge_snp_counts = get_snp_to_edge_counts( tree, snps_table_infile );
    set_tree_branch_length_to_snp_counts( edge_snp_counts, tree );
    // CommonTreeNewickWriter().write( tree, utils::to_file( tree_table_infile + "-brlen.newick" ));

    // Make a pairwise distance table for the tree.
    // make_node_distances_tables( tree, out_dir );

    // Get the sample data and compute per edge values for it
    auto const edge_values = make_sample_edge_snp_counts_derived(
        tree, sample_table_infile, exclude_damage
    );
    auto const clade_values = propagate_edge_snp_counts( tree, edge_values );

    // Write the table and the tree
    write_edge_snp_values( tree, edge_values, clade_values, out_base + ".csv" );
    sample_counts_tree_to_svg( tree, clade_values, out_base );
}

// =================================================================================================
//     Run Sample Dir
// =================================================================================================

std::vector<double> propagate_edge_values_inwards(
    Tree const& tree,
    std::vector<double> const& edge_values
) {
    LOG_DBG << "propagate_edge_values_inwards";
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

void make_sample_jplace(
    Tree const& tree,
    std::vector<double> const& edge_masses,
    std::string const& sample_name,
    std::string const& out_file
) {
    LOG_MSG << "make_sample_jplace";

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

void run_sample_dir(
    std::string const& tree_table_infile,
    std::string const& snps_table_infile,
    std::string const& sample_dir
) {
    // Get the output base name
    auto const out_dir = file_path( real_path( tree_table_infile )) + "/out/";

    // Get the tree
    auto tree = read_tree_from_table_file( tree_table_infile );

    // Set the branch lengths to the numbers of snps on each of them
    auto const edge_snp_counts = get_snp_to_edge_counts( tree, snps_table_infile );
    set_tree_branch_length_to_snp_counts( edge_snp_counts, tree );

    // We want an overview tree visualization of how many samples have their max on each branch.
    // Also, we do one where all clade values of the samples are accumulated, to see what happens.
    auto max_value_edge_counts = std::vector<double>( tree.edge_count(), 0 );
    auto accumulated_edge_values = std::vector<double>( tree.edge_count(), 0.0 );
    auto accumulated_clade_values = std::vector<double>( tree.edge_count(), 0.0 );

    // Process all sample files
    auto sample_files = dir_list_files( sample_dir, true, ".*calls$" );
    LOG_MSG << "found " << sample_files.size() << " sample files";
    for( size_t i = 0; i < sample_files.size(); ++i ) {
        LOG_MSG << "at sample " << i << ": " << sample_files[i];
        auto const out_base = out_dir + "/samples/" + file_basename( sample_files[i] );

        // Get the sample data and compute per edge values for it
        auto const edge_values = make_sample_edge_snp_counts_derived(
            tree, sample_files[i], false
        );
        auto const clade_values = propagate_edge_snp_counts( tree, edge_values );

        // Write the table and the tree.
        // Deactivated for now, as we have already produced those files.
        // write_edge_snp_values( tree, edge_values, clade_values, out_base + ".csv" );
        // sample_counts_tree_to_svg( tree, clade_values, out_base );

        // Make a fake jplace file from the sample
        make_sample_jplace(
            tree, clade_values, file_basename( sample_files[i] ), out_base + ".jplace"
        );

        // Accumulate the three vectors for the tree
        auto max_index = std::distance(
            clade_values.begin(), std::max_element( clade_values.begin(), clade_values.end() )
        );
        max_value_edge_counts[ max_index ] += 1.0;
        for( size_t j = 0; j < tree.edge_count(); ++j ) {
            accumulated_edge_values[j]  += edge_values[j];
            accumulated_clade_values[j] += clade_values[j];
        }
    }
    LOG_MSG << "done with samples";

    // Have done the below already, can deactivate for now.
    return;

    // Write values to lists for further examination.
    file_write( join( max_value_edge_counts, "\n" ) , out_dir + "max_value_edge_counts.txt" );
    file_write( join( accumulated_edge_values, "\n" ) , out_dir + "accumulated_edge_values.txt" );
    file_write( join( accumulated_clade_values, "\n" ) , out_dir + "accumulated_clade_values.txt" );

    // Write tree files of these values as well.
    sample_counts_tree_to_svg(
        tree, max_value_edge_counts, out_dir + "max_value_edge_counts_phylogram",
        CountSvgSettings( false, true, false )
    );
    sample_counts_tree_to_svg(
        tree, max_value_edge_counts, out_dir + "max_value_edge_counts_cladogram",
        CountSvgSettings( false, false, false )
    );
    sample_counts_tree_to_svg(
        tree, accumulated_clade_values, out_dir + "accumulated_clade_values",
        CountSvgSettings( false, true, false )
    );
    sample_counts_tree_to_svg(
        tree, accumulated_edge_values, out_dir + "accumulated_edge_values",
        CountSvgSettings( true, true, false )
    );

    // Accumulate values inwards
    auto const max_value_propagated_inwards = propagate_edge_values_inwards(
        tree, max_value_edge_counts
    );
    sample_counts_tree_to_svg(
        tree, max_value_propagated_inwards, out_dir + "max_value_propagated_inwards_phylogram",
        CountSvgSettings( false, true, false )
    );
    sample_counts_tree_to_svg(
        tree, max_value_propagated_inwards, out_dir + "max_value_propagated_inwards_cladogram",
        CountSvgSettings( false, false, false )
    );

    // Accumulate values outwards
    auto const max_value_propagated_outwards = propagate_edge_snp_counts(
        tree, max_value_edge_counts
    );
    sample_counts_tree_to_svg(
        tree, max_value_propagated_outwards, out_dir + "max_value_propagated_outwards_phylogram",
        CountSvgSettings( false, true, false )
    );
    sample_counts_tree_to_svg(
        tree, max_value_propagated_outwards, out_dir + "max_value_propagated_outwards_cladogram",
        CountSvgSettings( false, false, false )
    );

    // Also imbalances? Not really useful though.
    auto const max_value_imbalances = placement::epca_imbalance_vector( tree, max_value_edge_counts );
    sample_counts_tree_to_svg(
        tree, max_value_imbalances, out_dir + "max_value_imbalances",
        CountSvgSettings( true, true, false )
    );
}

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char const** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    Options::get().allow_file_overwriting( true );
    LOG_INFO << "started";

    // Get input
    if( argc != 4 ) {
        throw std::runtime_error( "Wrong usage" );
    }
    auto const tree_table_infile = std::string( argv[1] );
    auto const snps_table_infile = std::string( argv[2] );
    auto const sample_data       = std::string( argv[3] );

    // If the sample data is a single file, we run a vis of the sample.
    // If it's a directory, we instead run the process for all.
    if( is_file( sample_data )) {
        // Run with all or without damage
        run_single_sample( tree_table_infile, snps_table_infile, sample_data, false ) ;
        // run_single_sample( tree_table_infile, snps_table_infile, sample_data, true ) ;
    } else if( is_dir( sample_data )) {
        run_sample_dir( tree_table_infile, snps_table_infile, sample_data );
    } else {
        throw std::runtime_error( "Wrong usage" );
    }

    LOG_INFO << "finished";
    return 0;
}
