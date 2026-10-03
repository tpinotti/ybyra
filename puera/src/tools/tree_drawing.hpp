#ifndef PUERA_TOOLS_TREE_DRAWING_H_
#define PUERA_TOOLS_TREE_DRAWING_H_

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

#include "options/file_output.hpp"
#include "options/snps_table.hpp"
#include "options/tree_annotation.hpp"
#include "options/tree_output_svg.hpp"
#include "options/tree_table.hpp"

#include "genesis/tree/drawing/function.hpp"
#include "genesis/tree/tree.hpp"
#include "genesis/util/color/map.hpp"
#include "genesis/util/color/normalization.hpp"
#include "genesis/util/format/svg/group.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Tree Drawing
// =================================================================================================

/**
 * @brief Everything needed to draw a tree that is the same for all plots of a command.
 */
struct TreeDrawingSetup
{
    genesis::tree::Tree const&                   tree;
    genesis::tree::LayoutParameters              layout_params;
    double                                       size_unit;
    std::vector<genesis::util::format::SvgGroup> label_shapes;
};

/**
 * @brief Prepare the tree drawing, using the number of SNPs per branch as branch lengths,
 * if a SNP table was provided.
 */
TreeDrawingSetup prepare_tree_drawing(
    TreeTableOptions const&      tree_table,
    SnpsTableOptions const&      snps_table,
    SvgTreeOutputOptions const&  svg_tree,
    TreeAnnotationOptions const& annotation
);

/**
 * @brief Get the index of the edge leading to the placement @p node of a sample.
 */
size_t placement_edge_index(
    TreeTableOptions const& tree_table,
    std::string const&      node,
    std::string const&      sample_name
);

/**
 * @brief Draw the tree with branches colored by @p values, add the title, and write it
 * to the output file with the given @p infix.
 */
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
);

#endif // include guard
