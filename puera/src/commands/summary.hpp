#ifndef PUERA_COMMANDS_SUMMARY_H_
#define PUERA_COMMANDS_SUMMARY_H_

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

#include "options/file_output.hpp"

// =================================================================================================
//      Options
// =================================================================================================

class SummaryOptions
{
public:

    // Output options
    FileOutputOptions file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_summary( CLI::App& app );
void run_summary( SummaryOptions const& options );

#endif // include guard
