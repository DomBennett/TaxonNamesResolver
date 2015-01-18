#! /usr/bin/env python
# D.J. Bennett
# 01/06/2014
"""
TaxonNamesResolver is a python package for resolving taxonomic
names through Global Names Resolver (Copyright (C) 2012-2013
Marine Biological Laboratory).

Copyright (C) 2014  Dominic John Bennett

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
# Create namespace
from taxon_names_resolver.resolver import Resolver
from taxon_names_resolver.manip_tools import TaxDict
from taxon_names_resolver.manip_tools import taxTree
