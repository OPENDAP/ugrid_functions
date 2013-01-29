// -*- mode: c++; c-basic-offset:4 -*-

// This file is part of libdap, A C++ implementation of the OPeNDAP Data
// Access Protocol.

// Copyright (c) 2002,2003,2011,2012 OPeNDAP, Inc.
// Authors: Nathan Potter <ndp@opendap.org>
//          James Gallagher <jgallagher@opendap.org>
//          Scott Moe <smeest1@gmail.com>
//          Bill Howe <billhowe@cs.washington.edu>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// You can contact OPeNDAP, Inc. at PO Box 112, Saunderstown, RI. 02874-0112.

#ifndef _UgridUtilities_h
#define _UgridUtilities_h 1

#include <gridfields/array.h>


using namespace std;
using namespace libdap;

namespace ugrid {

class Array;

/**
 *  UGrid attribute vocabulary
 */
#define CF_ROLE "cf_role"
#define CF_STANDARD_NAME "standard_name"
#define UGRID_MESH_TOPOLOGY "mesh_topology"
#define UGRID_NODE_COORDINATES "node_coordinates"
#define UGRID_FACE_NODE_CONNECTIVITY "face_node_connectivity"
#define UGRID_DIMENSION "dimension"
#define UGRID_LOCATION "location"
#define UGRID_GRID_LOCATION "grid_location"
#define UGRID_NODE "node"
#define UGRID_MESH "mesh"
#define UGRID_START_INDEX "start_index"


GF::Array *extractGridFieldArray(libdap::Array *a, vector<int*> *sharedIntArrays, vector<float*> *sharedFloatArrays);
template<typename T> T *extract_array(libdap::Array * a);
template<typename DODS, typename T> T *extract_array_helper(libdap::Array *a);


string getAttributeValue(libdap::BaseType *bt, string aName) ;
bool matchesCfRoleOrStandardName(libdap::BaseType *bt, string aValue);
bool same_dimensions(libdap::Array *arr1, libdap::Array *arr2);

bool checkAttributeValue(libdap::BaseType *bt, string aName, string aValue);


vector<string> split(const string &s, char delim);
vector<string> &split(const string &s, char delim, vector<string> &elems);


int getNfrom3byNArray(libdap::Array *array);





}// namespace ugrid

#endif // _UgridUtilities_h
