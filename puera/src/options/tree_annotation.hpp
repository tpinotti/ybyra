#ifndef PUERA_OPTIONS_TREE_ANNOTATION_H_
#define PUERA_OPTIONS_TREE_ANNOTATION_H_

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

#include "options/tree_table.hpp"
#include "tools/cli_option.hpp"

#include "genesis/util/color/color.hpp"
#include "genesis/util/format/svg/document.hpp"
#include "genesis/util/format/svg/group.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Tree Annotation Options
// =================================================================================================

/**
 * @brief Options for the annotations drawn on top of a tree: clade labels, the placement markers,
 * and the title.
 *
 * All sizes are multipliers of automatic sizes that are relative to the size unit of the drawing,
 * see SvgTreeOutputOptions::size_unit().
 */
class TreeAnnotationOptions
{
public:

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_clade_labels_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Annotation"
    );

    void add_marker_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Annotation"
    );

    void add_title_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Annotation"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the clade labels as shapes per edge, to be drawn in the middle of the edge
     * leading to each labelled node. Empty if no labels file was provided.
     */
    std::vector<genesis::util::format::SvgGroup> make_clade_label_edge_shapes(
        TreeTableOptions const& tree_opts,
        double size_unit
    ) const;

    /**
     * @brief Get the marker ring to be drawn at a placement node. If @p dashed is set,
     * the ring is dashed, which we use to indicate that the placement resolves a score tie.
     */
    genesis::util::format::SvgGroup make_placement_marker(
        double size_unit,
        bool dashed
    ) const;

    /**
     * @brief Get a filled circle to be drawn at a node with placed samples. Its area is
     * proportional to @p fraction, which is the count at the node relative to the maximum count.
     */
    genesis::util::format::SvgGroup make_count_circle(
        double size_unit,
        double fraction
    ) const;

    /**
     * @brief Add a title line above the drawing, unless the title is deactivated.
     */
    void add_title(
        genesis::util::format::SvgDocument& doc,
        std::string const& title,
        double size_unit
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

private:

    genesis::util::color::Color marker_color_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> clade_labels_file_ = "";
    CliOption<double>      label_size_        = 1.0;

    CliOption<double>      marker_size_         = 1.0;
    CliOption<double>      marker_stroke_width_ = 1.0;
    CliOption<std::string> marker_color_opt_    = "#FF6600";

    CliOption<bool>        no_title_ = false;

};

#endif // include guard
