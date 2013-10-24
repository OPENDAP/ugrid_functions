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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//
// You can contact OPeNDAP, Inc. at PO Box 112, Saunderstown, RI. 02874-0112.

#ifndef _UgridUtilities_h
#define _UgridUtilities_h 1

#include <gridfields/array.h>


using namespace std;
using namespace libdap;

class libdap::Array;

namespace ugrid {

/**
 *  REQUIRED UGrid attribute vocabulary
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
#define UGRID_EDGE "edge"
#define UGRID_FACE "face"
#define UGRID_MESH "mesh"
#define UGRID_START_INDEX "start_index"

/**
 *  OPTIONAL UGrid attribute vocabulary
 */
#define UGRID_EDGE_NODE_CONNECTIVITY "edge_node_connectivity"


#define UGRID_FACE_COORDINATES "face_coordinates"
#define UGRID_EDGE_COORDINATES "edge_coordinates"
#define UGRID_FACE_EDGE_CONNECTIVITY "face_edge_connectivity"
#define UGRID_FACE_FACE_CONNECTIVITY "face_face_connectivity"

GF::Array *extractGridFieldArray(libdap::Array *a, vector<int*> *sharedIntArrays, vector<float*> *sharedFloatArrays);
GF::Array *newGFIndexArray(string name, long size, vector<int*> *sharedIntArrays);

/**
 * DAP Array data extraction helper method.
 */
template<typename DODS, typename T>T *extract_array_helper(libdap::Array *a) {
    int length = a->length();

    DODS *src = new DODS[length];

    a->value(src);

    T *dest = new T[length];

    for (int i = 0; i < length; ++i)
        dest[i] = (T) src[i];

    delete [] src;

    return dest;
}

/** Given a pointer to an Array that holds a numeric type, extract the
 values and return in an array of T. This function allocates the
 array using 'new T[n]' so delete[] can be used when you are done
 the data. */
template<typename T> T *extractArray(libdap::Array *a) {

    // Simple types are Byte, ..., Float64, String and Url.
    if ((a->type() == dods_array_c && !a->var()->is_simple_type())
            || a->var()->type() == dods_str_c || a->var()->type() == dods_url_c)
        throw Error(malformed_expr,
                "The function requires a DAP numeric-type array argument.");

    a->set_send_p(true);
    a->read();
    // This test should never pass due to the previous two lines;
    // reading here seems to make
    // sense rather than letting the caller forget to do so.
    // is read() idemopotent?
    if (!a->read_p())
        throw InternalErr(__FILE__, __LINE__,
                string("The Array '") + a->name()
                        + "'does not contain values. send_read_p() not called?");



    // The types of arguments that the CE Parser will build for numeric
    // constants are limited to Uint32, Int32 and Float64. See ce_expr.y.
    // Expanded to work for any numeric type so it can be used for more than
    // just arguments.
    switch (a->var()->type()) {
    case dods_byte_c:
        return extract_array_helper<dods_byte, T>(a);

    case dods_uint16_c:
        return extract_array_helper<dods_uint16, T>(a);

    case dods_int16_c:
        return extract_array_helper<dods_int16, T>(a);

    case dods_uint32_c:
        return extract_array_helper<dods_uint32, T>(a);

    case dods_int32_c:
        return extract_array_helper<dods_int32, T>(a);

    case dods_float32_c:
        // Added the following line. jhrg 8/7/12
        return extract_array_helper<dods_float32, T>(a);

    case dods_float64_c:
        return extract_array_helper<dods_float64, T>(a);

    default:
        throw InternalErr(__FILE__, __LINE__,
                "The argument list built by the CE parser contained an unsupported numeric type.");
    }
}

//template<typename T> T *extractArray(libdap::Array * a);
//template<typename DODS, typename T> T *extract_array_helper(libdap::Array *a);


string getAttributeValue(libdap::BaseType *bt, string aName) ;
bool matchesCfRoleOrStandardName(libdap::BaseType *bt, string aValue);

bool checkAttributeValue(libdap::BaseType *bt, string aName, string aValue);


vector<string> split(const string &s, char delim);
vector<string> &split(const string &s, char delim, vector<string> &elems);


int getNfrom3byNArray(libdap::Array *array);

libdap::Type getGridfieldsReturnType(libdap::Type type);




}// namespace ugrid

#endif // _UgridUtilities_h
