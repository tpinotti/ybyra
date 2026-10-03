#ifndef PUERA_OPTIONS_YBYRA_INPUT_H_
#define PUERA_OPTIONS_YBYRA_INPUT_H_

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

#include "options/file_input.hpp"
#include "options/tree_table.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <unordered_map>
#include <vector>

// =================================================================================================
//      Ybyra Input Options
// =================================================================================================

/**
 * @brief Options to read the per-sample results of a ybyra run.
 *
 * The sample calls are either given as individual files and directories, or as the ybyra output
 * directory, in which case the calls and placements are taken from their default locations in it.
 * Commands that only need the placements can leave out the calls options.
 */
class YbyraInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs
    // -------------------------------------------------------------------------

    struct SampleFile
    {
        std::string name;
        std::string path;
    };

    struct Placement
    {
        std::string node;
        double      tree_score;
        std::string flag;
    };

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add the input options. Without @p with_calls, only the placements are used,
     * which then are required.
     */
    void add_ybyra_input_opts_to_app(
        CLI::App* sub,
        bool with_calls = true,
        std::string const& group = "Input"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the sample names and paths of their calls files, sorted by name.
     *
     * The sample name is the file name without the `.calls` extension,
     * which is the same as the `individual` name used by ybyra in its placement tables.
     */
    std::vector<SampleFile> const& sample_files() const;

    /**
     * @brief Get the placements of the samples, by sample name.
     *
     * Empty if no placements file was provided, or if the ybyra directory does not contain one.
     */
    std::unordered_map<std::string, Placement> const& placements() const;

    /**
     * @brief Get the names of the samples that ybyra could not place, from the `fail.yplace`
     * file in the ybyra directory. Empty if not available.
     */
    std::vector<std::string> failed_samples() const;

    /**
     * @brief Read the calls file of a sample, and return the value per edge of the tree.
     *
     * Each derived call on an edge adds one, each ancestral call subtracts one.
     */
    std::vector<double> read_sample_edge_values(
        SampleFile const& sample,
        TreeTableOptions const& tree_opts
    ) const;

    /**
     * @brief Hint to give to users when the scores computed here do not match the ones from ybyra.
     */
    std::string mismatch_hint() const;

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

private:

    std::string placements_file_path_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    FileInputOptions       calls_input_;
    CliOption<std::string> ybyra_dir_;
    CliOption<std::string> placements_file_;
    CliOption<bool>        exclude_damage_ = false;
    bool                   with_calls_      = true;

    mutable std::vector<SampleFile> sample_files_;
    mutable std::unordered_map<std::string, Placement> placements_;
    mutable bool placements_read_ = false;

};

#endif // include guard
