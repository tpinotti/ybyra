#ifndef PUERA_OPTIONS_SNPS_TABLE_H_
#define PUERA_OPTIONS_SNPS_TABLE_H_

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

#include "options/tree_table.hpp"
#include "tools/cli_option.hpp"

#include "genesis/tree/tree.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// =================================================================================================
//      SNPs Table Options
// =================================================================================================

/**
 * @brief
 */
class SnpsTableOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SnpsTableOptions()  = default;
    virtual ~SnpsTableOptions() = default;

    SnpsTableOptions( SnpsTableOptions const& other ) = default;
    SnpsTableOptions( SnpsTableOptions&& )            = default;

    SnpsTableOptions& operator= ( SnpsTableOptions const& other ) = default;
    SnpsTableOptions& operator= ( SnpsTableOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief
     */
    void add_snps_table_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "SNPs Table"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Read the SNP table and return a vector with counts of SNPs for each edge.
     */
    std::vector<size_t> get_snps_per_edge_counts(
        TreeTableOptions const& tree_opts
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:


    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> snps_table_opt_ = "";

    CliOption<std::string> idx_col_opt_ = "id";
    CliOption<std::string> par_col_opt_ = "parent";
    CliOption<std::string> separator_char_opt_ = "tab";

};

#endif // include guard
