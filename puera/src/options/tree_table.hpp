#ifndef PUERA_OPTIONS_TREE_TABLE_H_
#define PUERA_OPTIONS_TREE_TABLE_H_

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

#include "tools/cli_option.hpp"

#include "genesis/tree/tree.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// =================================================================================================
//      Tree Options
// =================================================================================================

/**
 * @brief
 */
class TreeTableOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    TreeTableOptions()  = default;
    virtual ~TreeTableOptions() = default;

    TreeTableOptions( TreeTableOptions const& other ) = default;
    TreeTableOptions( TreeTableOptions&& )            = default;

    TreeTableOptions& operator= ( TreeTableOptions const& other ) = default;
    TreeTableOptions& operator= ( TreeTableOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief
     */
    void add_tree_table_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Tree Table"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the tree as provided by the user input table.
     */
    genesis::tree::Tree get_tree() const;

    /**
     * @brief Get a map of node names to edge indices in the given tree.
     */
    std::unordered_map<std::string, size_t> make_node_name_to_edge_index() const;

    /**
     * @brief Update the tree branch lengths to match the SNP count, for vis purposes.
     */
    void set_tree_branch_length_to_snp_counts(
        std::vector<size_t> const& edge_snp_counts
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:

    void read_tree_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> tree_table_opt_ = "";

    CliOption<std::string> idx_col_opt_ = "id";
    CliOption<std::string> par_col_opt_ = "parent";
    CliOption<std::string> separator_char_opt_ = "comma";

    mutable genesis::tree::Tree tree_;

};

#endif // include guard
