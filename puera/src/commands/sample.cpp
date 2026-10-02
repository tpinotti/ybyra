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

#include "commands/sample.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/util/core/logging.hpp"

#include <memory>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_sample( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SampleOptions>();
    auto sub = app.add_subcommand(
        "sample",
        "Plot the tree scores and placement of individual samples on the tree."
    );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( puera_cli_callback(
        sub,
        [ options ]() {
            run_sample( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_sample( SampleOptions const& options )
{
    (void) options;
    LOG_MSG << "The sample command is not yet implemented.";
}
