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

#include "options/tree_annotation.hpp"

#include "genesis/util/color/color.hpp"
#include "genesis/util/color/function.hpp"
#include "genesis/util/container/dataframe.hpp"
#include "genesis/util/container/dataframe/reader.hpp"
#include "genesis/util/format/svg/attribute.hpp"
#include "genesis/util/format/svg/shape.hpp"
#include "genesis/util/format/svg/text.hpp"
#include "genesis/util/io/input_source.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

using namespace genesis::util::color;
using namespace genesis::util::format;

// Sizes of the annotations, relative to the size unit of the drawing.
static double const label_font_factor_    = 0.037;
static double const marker_radius_factor_ = 0.03;
static double const marker_stroke_factor_ = 0.01;
static double const circle_radius_factor_ = 0.025;
static double const circle_stroke_factor_ = 0.003;
static double const circle_opacity_       = 0.7;

// =================================================================================================
//      Setup Functions
// =================================================================================================

void TreeAnnotationOptions::add_clade_labels_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    clade_labels_file_.option = sub->add_option(
        "--clade-labels-file",
        clade_labels_file_.value,
        "Tab-separated table with columns `node` and `label`, with a header line. "
        "Each listed node of the tree is labelled with its label text, "
        "drawn on the branch leading to the node. Typically used to label major haplogroups."
    );
    clade_labels_file_.option->check( CLI::ExistingFile );
    clade_labels_file_.option->group( group );

    label_size_.option = sub->add_option(
        "--label-size",
        label_size_.value,
        "Font size of the clade labels and the title, as a multiplier of the automatic size, "
        "which is scaled to the size of the tree."
    );
    label_size_.option->check( CLI::PositiveNumber );
    label_size_.option->group( group );
}

void TreeAnnotationOptions::add_marker_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    marker_size_.option = sub->add_option(
        "--marker-size",
        marker_size_.value,
        "Size of the markers that show sample placements, "
        "as a multiplier of the automatic size, which is scaled to the size of the tree."
    );
    marker_size_.option->check( CLI::PositiveNumber );
    marker_size_.option->group( group );

    marker_stroke_width_.option = sub->add_option(
        "--marker-stroke-width",
        marker_stroke_width_.value,
        "Stroke width of the markers that show sample placements, "
        "as a multiplier of the automatic width, which is scaled to the size of the tree."
    );
    marker_stroke_width_.option->check( CLI::PositiveNumber );
    marker_stroke_width_.option->group( group );

    marker_color_opt_.option = sub->add_option(
        "--marker-color",
        marker_color_opt_.value,
        "Color of the markers that show sample placements. "
        "Color can be specified in the format `#rrggbb` using hex values, or by web color names."
    );
    marker_color_opt_.option->group( group );
}

void TreeAnnotationOptions::add_title_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    no_title_.option = sub->add_flag(
        "--no-title",
        no_title_.value,
        "Do not add a title line above the tree."
    );
    no_title_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::vector<SvgGroup> TreeAnnotationOptions::make_clade_label_edge_shapes(
    TreeTableOptions const& tree_opts,
    double size_unit
) const {
    using namespace genesis::util::container;
    using namespace genesis::util::io;

    if( clade_labels_file_.value.empty() ) {
        return {};
    }
    auto const& path = clade_labels_file_.value;
    auto const& node_name_to_edge_index = tree_opts.node_name_to_edge_index();

    // Read the table and check its format.
    auto reader = DataframeReader<std::string>( '\t' ).row_names_from_first_col( false );
    auto const table = reader.read( from_file( path ));
    if( ! table.has_col_name( "node" ) || ! table.has_col_name( "label" )) {
        throw std::runtime_error(
            "Clade labels file " + path + " needs a header line with columns `node` and `label`"
        );
    }
    auto const& col_node  = table[ "node"  ].as<std::string>();
    auto const& col_label = table[ "label" ].as<std::string>();

    // Make a label for each listed node.
    auto const font_size = label_size_.value * label_font_factor_ * size_unit;
    auto result = std::vector<SvgGroup>( tree_opts.get_tree().edge_count() );
    std::unordered_set<std::string> seen;
    for( size_t i = 0; i < table.rows(); ++i ) {
        auto const& node  = col_node[i];
        auto const& label = col_label[i];
        if( seen.count( node ) > 0 ) {
            throw std::runtime_error(
                "Node \"" + node + "\" is listed multiple times in clade labels file " + path
            );
        }
        seen.insert( node );
        auto const edge_it = node_name_to_edge_index.find( node );
        if( edge_it == node_name_to_edge_index.end() ) {
            throw std::runtime_error(
                "Node \"" + node + "\" from clade labels file " + path + " is not in the tree, "
                "or is its root, which cannot be labelled. "
                "Is the labels file meant for this tree?"
            );
        }

        // Half-transparent white background circle, growing with the length of the label,
        // so that the label is readable on top of the branches.
        auto const radius = font_size * std::max( 0.5, 0.3 * label.size() + 0.1 );
        auto& group = result[ edge_it->second ];
        group.add( SvgCircle(
            SvgPoint( 0, 0 ), radius,
            SvgStroke( SvgStroke::Type::kNone ),
            SvgFill( Color( 1.0, 1.0, 1.0, 0.7 ))
        ));

        auto text = SvgText( label );
        text.font.size = font_size;
        text.anchor = SvgText::Anchor::kMiddle;
        text.dominant_baseline = SvgText::DominantBaseline::kCentral;
        text.dy = ".08em";
        group.add( text );
    }
    return result;
}

SvgGroup TreeAnnotationOptions::make_placement_marker(
    double size_unit,
    bool dashed
) const {
    auto stroke = SvgStroke( marker_color_(), marker_stroke_width_.value * marker_stroke_factor_ * size_unit );
    if( dashed ) {
        stroke.dash_array = { 1.5 * stroke.width, 1.0 * stroke.width };
    }

    SvgGroup result;
    result.add( SvgCircle(
        SvgPoint( 0, 0 ), marker_size_.value * marker_radius_factor_ * size_unit,
        stroke, SvgFill( SvgFill::Type::kNone )
    ));
    return result;
}

SvgGroup TreeAnnotationOptions::make_count_circle(
    double size_unit,
    double fraction
) const {
    auto color = marker_color_();
    color.a( circle_opacity_ );

    // White outline, so that overlapping circles stay distinguishable.
    auto const radius = marker_size_.value * circle_radius_factor_ * size_unit * std::sqrt( fraction );
    SvgGroup result;
    result.add( SvgCircle(
        SvgPoint( 0, 0 ), radius,
        SvgStroke( Color( 1.0, 1.0, 1.0 ), marker_stroke_width_.value * circle_stroke_factor_ * size_unit ),
        SvgFill( color )
    ));
    return result;
}

void TreeAnnotationOptions::add_title(
    SvgDocument& doc,
    std::string const& title,
    double size_unit
) const {
    if( no_title_.value ) {
        return;
    }

    // Place the title above the top left corner of the current drawing.
    auto const font_size = label_size_.value * label_font_factor_ * size_unit;
    auto const bbox = doc.bounding_box();
    auto text = SvgText( title, SvgPoint( bbox.top_left.x, bbox.top_left.y - font_size ));
    text.font.size = font_size;
    doc.add( text );
}

// =================================================================================================
//      Internal Helpers
// =================================================================================================

Color TreeAnnotationOptions::marker_color_() const
{
    try {
        return resolve_color_string( marker_color_opt_.value );
    } catch( std::exception const& ex ) {
        throw CLI::ValidationError(
            "--marker-color", "Invalid color '" + marker_color_opt_.value + "': " + ex.what()
        );
    }
}
