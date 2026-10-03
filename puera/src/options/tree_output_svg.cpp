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

#include "options/tree_output_svg.hpp"

#include <algorithm>

// Branch stroke width, relative to the size unit of the drawing.
static double const stroke_width_factor_ = 0.004;

// =================================================================================================
//      Setup Functions
// =================================================================================================

void SvgTreeOutputOptions::add_svg_tree_output_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Shape
    shape_.option = sub->add_option(
        "--svg-tree-shape",
        shape_.value,
        "Shape of the tree."
    );
    shape_.option->group( group );
    shape_.option->transform(
        CLI::IsMember({ "circular", "rectangular" }, CLI::ignore_case )
    );

    // Type
    type_.option = sub->add_option(
        "--svg-tree-type",
        type_.value,
        "Type of the tree, either using branch lengths (`phylogram`), or not (`cladogram`). "
        "With `auto`, a phylogram is drawn if a SNP table is provided, using the number of SNPs "
        "on each branch as its length, and a cladogram otherwise."
    );
    type_.option->group( group );
    type_.option->transform(
        CLI::IsMember({ "auto", "cladogram", "phylogram" }, CLI::ignore_case )
    );

    // Stroke width
    stroke_width_.option = sub->add_option(
        "--svg-tree-stroke-width",
        stroke_width_.value,
        "Stroke width of the branches of the tree, as a multiplier of the automatic width, "
        "which is scaled to the size of the tree."
    );
    stroke_width_.option->group( group );
    stroke_width_.option->check( CLI::PositiveNumber );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::tree::LayoutParameters SvgTreeOutputOptions::layout_parameters(
    genesis::tree::Tree const& tree,
    bool has_branch_lengths
) const {
    using namespace genesis::tree;
    LayoutParameters result;

    if( shape_.value == "circular" ) {
        result.shape = LayoutShape::kCircular;
    } else if( shape_.value == "rectangular" ) {
        result.shape = LayoutShape::kRectangular;
    } else {
        throw CLI::ValidationError( "--svg-tree-shape", "Invalid shape '" + shape_.value + "'." );
    }

    if( type_.value == "phylogram" && ! has_branch_lengths ) {
        throw CLI::ValidationError(
            "--svg-tree-type", "Drawing a phylogram requires branch lengths, "
            "which are computed from the SNP table. Please provide a SNP table."
        );
    }
    if( type_.value == "phylogram" || ( type_.value == "auto" && has_branch_lengths )) {
        result.type = LayoutType::kPhylogram;
    } else if( type_.value == "cladogram" || type_.value == "auto" ) {
        result.type = LayoutType::kCladogram;
    } else {
        throw CLI::ValidationError( "--svg-tree-type", "Invalid type '" + type_.value + "'." );
    }

    result.stroke.width = stroke_width_.value * stroke_width_factor_ * size_unit( tree );
    result.ladderize = true;
    return result;
}

double SvgTreeOutputOptions::size_unit( genesis::tree::Tree const& tree )
{
    // Same as the automatic radius of the genesis CircularLayout.
    return std::max( 50.0, static_cast<double>( tree.node_count() ));
}
