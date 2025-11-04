#ifndef PUERA_TOOLS_CLI_SETUP_H_
#define PUERA_TOOLS_CLI_SETUP_H_

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

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

// =================================================================================================
//      CLI11 Setup
// =================================================================================================

/**
 * @brief Callback wrapper function to use for puera commands.
 *
 * Using this function to generate the callback for CLI commands ensures a consistent interface:
 * It takes the subcommand, so that its options can be printed, as well as the actual callback
 * function, identically to what CLI expects.
 */
std::function<void()> puera_cli_callback(
    CLI::App const*          subcommand,
    std::function<void()>    run_function
);

// =================================================================================================
//      Checks and Helpers
// =================================================================================================

/**
 * @brief Check recursively that all puera subcommands have unique names, across modules.
 *
 * We want to ensure clarity of the commands, so we do not allow identically named subcommands
 * in different modules.
 */
void check_unique_command_names( CLI::App const& app );

/**
 * @brief Check recursively that all puera subcommands and options have names and descriptions set.
 */
void check_subcommand_names( CLI::App const& app );

/**
 * @brief Fix/capture the current default values of all options of an app, recursively.
 *
 * In release 1.8 of CLI, the default capturing mechanism for variables was changed. Previously,
 * it was done via a `true` flag in the `add_option` function, but now, `capture_default_str()`
 * is the new way. In order to make sure that this is called for every option of every command,
 * run this function.
 *
 * This serves also a second purpose: Once the current values of the variables have been captured,
 * they are availalbe in the help messages, and (more importantly) also can be used for our output
 * of current option settings that we print when running a command.
 * See puera_cli_callback() for details.
 */
void fix_cli_default_values( CLI::App& app );

/**
 * @brief Set the help group name for all subcommands of a module.
 */
void set_module_help_group( CLI::App& module, std::string const& group_name = "Global Options" );

#endif // include guard
