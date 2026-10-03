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

#include "options/color_norm.hpp"

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* ColorNormOptions::add_scaling_opt_to_app(
    CLI::App* sub,
    std::string const& auto_description,
    std::string const& group
) {
    bool const with_auto = ! auto_description.empty();
    if( with_auto ) {
        scaling_option.value = "auto";
    }
    scaling_option = sub->add_option(
        "--color-scaling",
        scaling_option.value,
        "Scaling of the color scale, `linear` or logarithmic (`log`). " + auto_description
    );
    scaling_option.option->group( group );
    if( with_auto ) {
        scaling_option.option->transform(
            CLI::IsMember({ "auto", "linear", "log" }, CLI::ignore_case )
        );
    } else {
        scaling_option.option->transform(
            CLI::IsMember({ "linear", "log" }, CLI::ignore_case )
        );
    }
    return scaling_option.option;
}

CLI::Option* ColorNormOptions::add_min_value_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    min_value_option = sub->add_option(
        "--min-value",
        min_value_option.value,
        "Minimum value that is represented by the color scale. "
        "If not set, a default suitable for the plot is used."
    );
    min_value_option.option->default_function( [](){ return std::string( "auto" ); } );
    min_value_option.option->group( group );
    return min_value_option.option;
}

CLI::Option* ColorNormOptions::add_max_value_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    max_value_option = sub->add_option(
        "--max-value",
        max_value_option.value,
        "Maximum value that is represented by the color scale. "
        "If not set, the maximum value of the data is used."
    );
    max_value_option.option->default_function( [](){ return std::string( "auto" ); } );
    max_value_option.option->group( group );
    return max_value_option.option;
}

CLI::Option* ColorNormOptions::add_mask_value_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    mask_value_option = sub->add_option(
        "--mask-value",
        mask_value_option.value,
        "Values of the data that compare equal to this mask value are colored using "
        "`--mask-color`. Infinities and NaN values are always masked. "
        "If not set, no mask value is applied."
    );
    mask_value_option.option->group( group );
    return mask_value_option.option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::unique_ptr<genesis::util::color::ColorNormalizationLinear>
ColorNormOptions::get_sequential_norm( bool auto_log ) const
{
    using namespace genesis::util::color;

    std::unique_ptr<ColorNormalizationLinear> result;
    if( log_scaling( auto_log )) {
        result = std::make_unique<ColorNormalizationLogarithmic>();
    } else {
        result = std::make_unique<ColorNormalizationLinear>();
    }
    apply_options( *result );
    return result;
}

void ColorNormOptions::apply_options( genesis::util::color::ColorNormalizationLinear& norm ) const
{
    if( min_value_option.is_set() ) {
        norm.min_value( min_value_option.value );
    }
    if( max_value_option.is_set() ) {
        norm.max_value( max_value_option.value );
    }
    if( mask_value_option.is_set() ) {
        norm.mask_value( mask_value_option.value );
    }
}
