#ifndef PUERA_OPTIONS_COLOR_NORM_H_
#define PUERA_OPTIONS_COLOR_NORM_H_

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

#include "genesis/util/color/norm_linear.hpp"
#include "genesis/util/color/norm_logarithmic.hpp"

#include <limits>
#include <memory>
#include <string>

// =================================================================================================
//      Color Norm Options
// =================================================================================================

/**
 * @brief Helper class to add command line options for a sequential color normalization.
 *
 * The values set by the user take precedence over the values that a command sets on the
 * normalization itself, see apply_options().
 */
class ColorNormOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    ColorNormOptions()  = default;
    ~ColorNormOptions() = default;

    ColorNormOptions( ColorNormOptions const& other ) = default;
    ColorNormOptions( ColorNormOptions&& )            = default;

    ColorNormOptions& operator= ( ColorNormOptions const& other ) = default;
    ColorNormOptions& operator= ( ColorNormOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add the option to choose linear or logarithmic scaling. If @p auto_description is
     * given, `auto` is also offered as the default, letting the command choose the scaling,
     * as explained by the description.
     */
    CLI::Option* add_scaling_opt_to_app(
        CLI::App* sub,
        std::string const& auto_description = "",
        std::string const& group = "Color"
    );

    CLI::Option* add_min_value_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Color"
    );

    CLI::Option* add_max_value_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Color"
    );

    CLI::Option* add_mask_value_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Color"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Return whether to use logarithmic scaling, with @p auto_log used for `auto`.
     */
    bool log_scaling( bool auto_log = false ) const
    {
        return scaling_option.value == "log" || ( scaling_option.value == "auto" && auto_log );
    }

    bool min_value_is_set() const
    {
        return min_value_option.is_set();
    }

    bool max_value_is_set() const
    {
        return max_value_option.is_set();
    }

    /**
     * @brief Get a linear or logarithmic normalization, depending on the log scaling option,
     * with the user provided values already applied.
     */
    std::unique_ptr<genesis::util::color::ColorNormalizationLinear> get_sequential_norm(
        bool auto_log = false
    ) const;

    /**
     * @brief Overwrite the min, max, and mask values of the @p norm with the values
     * provided by the user, if any.
     */
    void apply_options( genesis::util::color::ColorNormalizationLinear& norm ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> scaling_option = "linear";
    CliOption<double>      min_value_option  = 0.0;
    CliOption<double>      max_value_option  = 1.0;
    CliOption<double>      mask_value_option = std::numeric_limits<double>::quiet_NaN();

};

#endif // include guard
