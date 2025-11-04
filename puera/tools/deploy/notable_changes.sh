#!/bin/bash

# puera - visualizing y-chromosome samples on trees
# Copyright (C) 2025 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@sund.ku.dk>
# University of Copenhagen, Globe Institute, Section for GeoGenetics
# Oster Voldgade 5-7, 1350 Copenhagen K, Denmark

####################################################################################################
#    This script lists notable changes since the last release.
####################################################################################################

# Get current version.
last_tag=`git describe --abbrev=0 --tags`

# Get all important commits after the last tag and format them for Markdown.
echo -e "\e[34mNotable changes since version ${last_tag}\e[0m\n"
git log ${last_tag}..HEAD --oneline | cut -d " " -f 1 --complement | egrep -iv "^(Minor|Merge|Release)" | sed "s/^/  \* /g"
