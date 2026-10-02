#ifndef PUERA_OPTIONS_SAMPLES_H_
#define PUERA_OPTIONS_SAMPLES_H_

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

#include "options/file_input.hpp"
#include "options/tree_table.hpp"
#include "tools/cli_option.hpp"

#include "genesis/tree/tree.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// =================================================================================================
//      Tree Table Options
// =================================================================================================

/**
 * @brief
 */
class SamplesOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SamplesOptions()  = default;
    virtual ~SamplesOptions() = default;

    SamplesOptions( SamplesOptions const& other ) = default;
    SamplesOptions( SamplesOptions&& )            = default;

    SamplesOptions& operator= ( SamplesOptions const& other ) = default;
    SamplesOptions& operator= ( SamplesOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief
     */
    void add_samples_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Samples"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the tree as provided by the user input table.
     */
    std::vector<double> make_sample_edge_snp_counts_derived(
        TreeTableOptions const& tree_opts,
        size_t smp_idx
    );

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:


    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    FileInputOptions samples_opt_;

    // CliOption<std::string> snp_col_opt_ = "snpId";
    CliOption<std::string> idx_col_opt_ = "id";
    CliOption<std::string> par_col_opt_ = "parent";
    CliOption<std::string> stt_col_opt_ = "state";
    CliOption<std::string> dmg_col_opt_ = "damage";

    CliOption<std::string> separator_char_opt_ = "tab";

    CliOption<bool> exclude_damage_ = false;

};

#endif // include guard
