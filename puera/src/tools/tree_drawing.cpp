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

#include "tools/tree_drawing.hpp"

#include <stdexcept>

// =================================================================================================
//      Tree Drawing
// =================================================================================================

TreeDrawingSetup prepare_tree_drawing(
    TreeTableOptions const&      tree_table,
    SnpsTableOptions const&      snps_table,
    SvgTreeOutputOptions const&  svg_tree,
    TreeAnnotationOptions const& annotation
) {
    if( snps_table.provided() ) {
        tree_table.set_tree_branch_length_to_snp_counts(
            snps_table.get_snps_per_edge_counts( tree_table )
        );
    }
    auto const& tree = tree_table.get_tree();
    auto const size_unit = SvgTreeOutputOptions::size_unit( tree );
    return TreeDrawingSetup{
        tree,
        svg_tree.layout_parameters( tree, snps_table.provided() ),
        size_unit,
        annotation.make_clade_label_edge_shapes( tree_table, size_unit )
    };
}

size_t placement_edge_index(
    TreeTableOptions const& tree_table,
    std::string const&      node,
    std::string const&      sample_name
) {
    auto const& node_name_to_edge_index = tree_table.node_name_to_edge_index();
    auto const edge_it = node_name_to_edge_index.find( node );
    if( edge_it == node_name_to_edge_index.end() ) {
        throw std::runtime_error(
            "Placement node \"" + node + "\" of sample " + sample_name + " is not in the tree, "
            "or is its root. Is the tree table the same as the one used in the ybyra run?"
        );
    }
    return edge_it->second;
}

void write_tree_svg(
    TreeDrawingSetup const&                             setup,
    std::vector<double> const&                          values,
    genesis::util::color::ColorMap const&               color_map,
    genesis::util::color::ColorNormalization const&     color_norm,
    std::vector<genesis::util::format::SvgGroup> const& node_shapes,
    TreeAnnotationOptions const&                        annotation,
    std::string const&                                  title,
    FileOutputOptions const&                            file_output,
    std::string const&                                  infix
) {
    auto doc = genesis::tree::get_color_tree_svg_document(
        setup.tree, setup.layout_params, color_map( color_norm, values ), color_map, color_norm,
        node_shapes, setup.label_shapes
    );
    annotation.add_title( doc, title, setup.size_unit );
    doc.write( file_output.get_output_target( infix, "svg" ));
}
