#ifndef PUERA_FUNCTIONS_FUNCTIONS_H_
#define PUERA_FUNCTIONS_FUNCTIONS_H_

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

#include "genesis/tree/tree.hpp"

#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// =================================================================================================
//      Helper Functions
// =================================================================================================

std::vector<double> propagate_edge_values_inwards(
    genesis::tree::Tree const& tree,
    std::vector<double> const& edge_values
);

std::vector<double> propagate_edge_snp_counts(
    genesis::tree::Tree const& tree,
    std::vector<double> const& edge_snp_counts
);

#endif // include guard
