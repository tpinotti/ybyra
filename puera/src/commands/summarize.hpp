#ifndef PUERA_COMMANDS_SUMMARIZE_H_
#define PUERA_COMMANDS_SUMMARIZE_H_

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

#include "CLI/CLI.hpp"

#include "options/color_map.hpp"
#include "options/file_input.hpp"
#include "options/file_output.hpp"
#include "options/samples.hpp"
#include "options/snps_table.hpp"
#include "options/tree_table.hpp"
#include "tools/cli_option.hpp"

#include <limits>
#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class SummarizeOptions
{
public:

    // Input options
    TreeTableOptions tree_table_input;
    SnpsTableOptions snps_table_input;
    SamplesOptions   sample_input;

    // Plot settings
    ColorMapOptions color_map;
    CliOption<std::string> color_norm = "linear";
    // CliOption<double> min_value = std::numeric_limits<double>::quiet_NaN();
    // CliOption<double> max_value = std::numeric_limits<double>::quiet_NaN();

    // Output options
    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_summarize( CLI::App& app );
void run_summarize( SummarizeOptions const& options );

#endif // include guard
