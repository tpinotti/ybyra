#ifndef PUERA_OPTIONS_TREE_OUTPUT_SVG_H_
#define PUERA_OPTIONS_TREE_OUTPUT_SVG_H_

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

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"

#include "genesis/tree/drawing/function.hpp"
#include "genesis/tree/tree.hpp"

#include <string>

// =================================================================================================
//      SVG Tree Output Options
// =================================================================================================

/**
 * @brief Options for the shape, type, and stroke width of SVG tree drawings.
 */
class SvgTreeOutputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_svg_tree_output_opts_to_app(
        CLI::App* sub,
        std::string const& group = "SVG Tree"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the layout parameters for drawing the @p tree.
     *
     * With the `auto` tree type, a phylogram is drawn if @p has_branch_lengths is set,
     * and a cladogram otherwise.
     */
    genesis::tree::LayoutParameters layout_parameters(
        genesis::tree::Tree const& tree,
        bool has_branch_lengths
    ) const;

    /**
     * @brief Get the base size unit of the drawing, to which all other sizes are relative.
     *
     * Genesis uses the number of nodes in the tree as the radius of circular trees, and roughly
     * 2*pi times that as the height of rectangular trees. Line widths, font sizes etc that are
     * proportional to this unit hence look the same for trees of all sizes and both shapes.
     */
    static double size_unit( genesis::tree::Tree const& tree );

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> shape_        = "circular";
    CliOption<std::string> type_         = "auto";
    CliOption<double>      stroke_width_ = 1.0;

};

#endif // include guard
