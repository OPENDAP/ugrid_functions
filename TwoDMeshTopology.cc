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

#include "config.h"

#include <sstream>
#include <vector>

#include <gridfields/type.h>
#include <gridfields/gridfield.h>
#include <gridfields/grid.h>
#include <gridfields/implicit0cells.h>
#include <gridfields/gridfieldoperator.h>
#include <gridfields/restrict.h>
#include <gridfields/refrestrict.h>
#include <gridfields/array.h>


#include "BaseType.h"
#include "Int32.h"
#include "Float64.h"
#include "Array.h"
#include "util.h"


#include "ugrid_utils.h"
#include "NDimensionalArray.h"
#include "MeshDataVariable.h"
#include "TwoDMeshTopology.h"

#include "BESDebug.h"

using namespace std;
using std::endl;
using namespace libdap;
using namespace ugrid;

namespace ugrid {


TwoDMeshTopology::TwoDMeshTopology():myVar(0),gridTopology(0),d_inputGridField(0),resultGridField(0),sharedNodeArray(0)
{
    rangeDataArrays = new vector<MeshDataVariable *>();
    sharedIntArrays = new vector<int *>();
    sharedFloatArrays = new vector<float *>();

    nodeCount = 0;
    nodeCoordinateArrays = 0;
    faceNodeConnectivityArray = 0;
    faceCount = 0;
    faceCoordinateNames = 0;
    faceCoordinateArrays = 0;
    _initialized = false;
}


TwoDMeshTopology::~TwoDMeshTopology()
{
	BESDEBUG("ugrid", "~TwoDMeshTopology() - BEGIN" << endl);
	BESDEBUG("ugrid", "~TwoDMeshTopology() - (" << this << ")" << endl);

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting resultGridField." << endl);
	delete resultGridField;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting inputGridField." << endl);
	delete d_inputGridField;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting gridTopology." << endl);
	delete gridTopology;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting sharedIntArrays..." << endl);
	for (vector<int *>::iterator it = sharedIntArrays->begin(); it != sharedIntArrays->end(); ++it) {
		int *i = *it;
		BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting int array '" << i << "'" << endl);
		delete [] i;
	}
	delete sharedIntArrays;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting sharedFloatArrays..." << endl);
	for (vector<float *>::iterator it = sharedFloatArrays->begin(); it != sharedFloatArrays->end(); ++it) {
		float *f = *it;
		BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting float array '" << f << "'" << endl);
		delete [] f;
	}
	delete sharedFloatArrays;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting MeshDataVariables..." << endl);
	vector<MeshDataVariable *>::iterator mdvIt;
	for (mdvIt = rangeDataArrays->begin(); mdvIt != rangeDataArrays->end(); ++mdvIt) {
		MeshDataVariable *mdv = *mdvIt;
		BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting MeshDataVariable '"<< mdv->getName() << "'" << endl);
		delete mdv;
	}
	delete rangeDataArrays;

	BESDEBUG("ugrid", "~TwoDMeshTopology() - Deleting GF::Node array." << endl);
	delete sharedNodeArray;


	BESDEBUG("ugrid", "~TwoDMeshTopology() - END" << endl);
}


/**
 * @TODO only call this from the constructor?? Seems like the thing to do, but then the constructor will be possibly throwing an Error - Is that OK?
 */
void TwoDMeshTopology::init(string meshVarName, DDS &dds)
{

    if(_initialized)
        return;

	myVar = dds.var(meshVarName);

	dimension = getAttributeValue(myVar,UGRID_DIMENSION);

    if (dimension.empty()) {
        string msg = "TwoDMeshTopology::init() - The mesh topology variable  '" + myVar->name() + "' is missing the required attribute named '" + UGRID_DIMENSION + "'";
        BESDEBUG("ugrid",msg );
        throw Error( msg);
    }

	// Retrieve the node coordinate arrays for the mesh
	ingestNodeCoordinateArrays(myVar, dds);


    // Why would we retrieve this? I think Bill said that this needs to be recomputed after a restrict operation.
    //  @TODO Verify that Bill actually said that this needs to be recomputed.
    // Retrieve the face coordinate arrays (if any)  for the mesh
    ingestFaceCoordinateArrays(myVar,dds);


    // Inspect and QC the face node connectivity array for the mesh
    ingestFaceNodeConnectivityArray(myVar,dds);



	_initialized = true;

}

void TwoDMeshTopology::setNodeCoordinateDimension(MeshDataVariable *mdv)
{
    BESDEBUG("ugrid", "TwoDMeshTopology::setNodeCoordinateDimension() - BEGIN" << endl);
    libdap::Array *dapArray = mdv->getDapArray();
    libdap::Array::Dim_iter ait1;


    for (ait1 = dapArray->dim_begin();   ait1 != dapArray->dim_end();      ++ait1) {

        string dimName = dapArray->dimension_name(ait1);

        if(dimName.compare(nodeDimensionName) == 0){ // are the names the same?
            BESDEBUG("ugrid", "TwoDMeshTopology::setNodeCoordinateDimension() - Found dimension name matching nodeDimensionName '"<< nodeDimensionName << "'" << endl);
            int size = dapArray->dimension_size(ait1,true);
            if(size == nodeCount){ // are they the same size?
                BESDEBUG("ugrid", "TwoDMeshTopology::setNodeCoordinateDimension() - Dimension sizes match (" << libdap::long_to_string(nodeCount) << ") - DONE" << endl);
                // Yay! We found the node coordinate dimension
                mdv->setLocationCoordinateDimension(ait1);
                return;
            }
        }

    }

    throw Error(
            "Unable to determine the node coordinate dimension of the range variable '"
                    + mdv->getName() + "'  The node dimension is named '"+ nodeDimensionName
                    + "'  with size " + libdap::long_to_string(nodeCount)
    );


}


void TwoDMeshTopology::setFaceCoordinateDimension(MeshDataVariable *mdv)
{


    libdap::Array *dapArray = mdv->getDapArray();
    libdap::Array::Dim_iter ait1;

    for (ait1 = dapArray->dim_begin();   ait1 != dapArray->dim_end();      ++ait1) {
        string dimName = dapArray->dimension_name(ait1);

        if(dimName.compare(faceDimensionName) == 0){ // are the names the same?
            int size = dapArray->dimension_size(ait1,true);
            if(size == faceCount){ // are they the same size?
                // Yay! We found the node coordinate dimension
                mdv->setLocationCoordinateDimension(ait1);
                return;
            }
        }

    }
    throw Error(
            "Unable to determine the face coordinate dimension of the range variable '"
                    + mdv->getName() + "'  The face coordinate dimension is named '"+ faceDimensionName
                    + "' with size " + libdap::long_to_string(faceCount)
    );



}


void TwoDMeshTopology::addDataVariable(MeshDataVariable *mdv)
{

    BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - BEGIN" << endl);

    libdap::Array *dapArray = mdv->getDapArray();

    switch(mdv->getGridLocation()){

    case node:
    {
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Checking Node variable  '"<< mdv->getName() << "'" << endl);

        // Locate and set the MDV's node coordinate dimension.
        setNodeCoordinateDimension(mdv);
        rangeDataArrays->push_back(mdv);
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);
    }
    break;

    case edge:
    {
        //@TODO Implement sanity checks... some kind of input QC.
        rangeDataArrays->push_back(mdv);
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);
    }
    break;

    case face:
    {
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Checking Face variable  '"<< mdv->getName() << "'" << endl);

        // Locate and set the MDV's face coordinate dimension.
        setFaceCoordinateDimension(mdv);
        rangeDataArrays->push_back(mdv);
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);

    }
    break;

    default:
    {
        string msg = "TwoDMeshTopology::addDataVariable() - Unknown location value '"  + libdap::long_to_string(mdv->getGridLocation()) + "'";
        BESDEBUG("ugrid",  msg << endl);
        throw Error( msg );
    }
    break;
    }

    BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - DONE" << endl);

#if 0
    BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - BEGIN" << endl);

	libdap::Array *dapArray = mdv->getDapArray();

	switch(mdv->getGridLocation()){

	case node:
	{
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Checking Node variable  '"<< mdv->getName() << "'" << endl);
	    // Make sure that the requested range variable is the same shape as the node coordinate arrays
	    // We only need to test the first nodeCoordinate array against the first rangeVar array
	    // because we have already made sure all of the node coordinate arrays are the same size and
	    // that all the rangeVar arrays are the same size. This is just to compare the two collections.
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Comparing dimensionality of '"<< mdv->getName() << "' to node coordinate arrays." << endl);
	    libdap::Array *firstCoordinate = (*(nodeCoordinateArrays))[0];
	     if (!same_dimensions(dapArray, firstCoordinate)){
	         throw Error(
	                 "The dimensions of the requested range variable "
	                         + mdv->getName() + " has a location of 'node' (aka rank 0) but it does not match the shape "
	                         + " of the node coordinate array "
	                         + firstCoordinate->name());
	     }
	     rangeDataArrays->push_back(mdv);
         BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);
	}
	break;

	case edge:
	{
	    //@TODO Implement sanity checks... some kind of input QC.
        rangeDataArrays->push_back(mdv);
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);
	}
	break;

	case face:
	{
        BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Checking Face variable  '"<< mdv->getName() << "'" << endl);


	    // Check face coordinate arrays, if present, to be sure the sizes match.
	    if(faceCoordinateArrays!=0 && !faceCoordinateArrays->empty()){
	         BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Comparing dimensionality of '"<< mdv->getName() << "' to face coordinate arrays." << endl);
            libdap::Array *firstCoordinate = (*(faceCoordinateArrays))[0];
             if (!same_dimensions(dapArray, firstCoordinate)){
                 throw Error(
                         "The dimensions of the requested range variable "
                                 + mdv->getName() + " has a location of 'face' (aka rank of 2) but it does not match the shape "
                                 + " of the face coordinate array "
                                 + firstCoordinate->name());
             }
	     }

         rangeDataArrays->push_back(mdv);
         BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - Added MeshDataVariable '"<< mdv->getName() << "' " << endl);

	}
	break;

	default:
	{
	    string msg = "TwoDMeshTopology::addDataVariable() - Unknown location value '"  + libdap::long_to_string(mdv->getGridLocation()) + "'";
	    BESDEBUG("ugrid",  msg << endl);
        throw Error( msg );
	}
    break;
	}

    BESDEBUG("ugrid", "TwoDMeshTopology::addDataVariable() - DONE" << endl);
#endif


}




/**
 * Locates the the DAP variable identified by the face_node_connectivity attribute of the
 * meshTopology variable.
 */
void TwoDMeshTopology::ingestFaceNodeConnectivityArray(libdap::BaseType *meshTopology, libdap::DDS &dds)
{

	BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Locating FNCA" << endl);

	string face_node_connectivity_var_name;
    AttrTable at = meshTopology->get_attr_table();

    AttrTable::Attr_iter iter_fnc = at.simple_find(UGRID_FACE_NODE_CONNECTIVITY);
    if (iter_fnc != at.attr_end()) {
    	face_node_connectivity_var_name = at.get_attr(iter_fnc, 0);
    }
    else {
    	throw Error("Could not locate the "  UGRID_FACE_NODE_CONNECTIVITY  " attribute in the "   UGRID_MESH_TOPOLOGY   " variable! "
    			"The mesh_topology variable is named " + meshTopology->name());
    }
    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Located the '" << UGRID_FACE_NODE_CONNECTIVITY << "' attribute." << endl);

	// Find the variable using the name

    BaseType *btp = dds.var(face_node_connectivity_var_name);

    if(btp==0)
    	throw Error("Could not locate the " UGRID_FACE_NODE_CONNECTIVITY " variable named '" + face_node_connectivity_var_name + "'! "+
    			"The mesh_topology variable is named "+meshTopology->name());


    // Is it an array?
    libdap::Array *fncArray = dynamic_cast<libdap::Array*>(btp);
    if(fncArray == 0) {
        throw Error(malformed_expr,"Face Node Connectivity variable '"+face_node_connectivity_var_name+"' is not an Array type. It's an instance of " + btp->type_name());
    }

    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Located FNC Array '" << fncArray->name() << "'." << endl);

    // It's got to have exactly 2 dimensions - [max#nodes_per_face][#faces]
    int numDims = fncArray->dimensions(true) ;
    if( numDims != 2){
        throw Error(malformed_expr,"Face Node Connectivity variable '"+face_node_connectivity_var_name
                 +"' Must have two (2) dimensions. It has " + libdap::long_to_string(numDims));
     }

    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - FNC Array '" << fncArray->name() << "' has two (2) dimensions." << endl);


    // We just need to have both dimensions handy, so get the first and the second
    libdap::Array::Dim_iter firstDim = fncArray->dim_begin();
    libdap::Array::Dim_iter secondDim = fncArray->dim_begin(); // same as the first for a moment...
    secondDim++; // now it's second!


    if(faceDimensionName.empty()){
        // By now we know it only has two dimensions, but since there is no promise that they'll be in a particular order
        // we punt: We'll assume that smallest of the two is in fact the nodes per face and the larger the face index dimensions.
        int sizeFirst  = fncArray->dimension_size(firstDim,true);
        int sizeSecond = fncArray->dimension_size(secondDim,true);
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - sizeFirst: "<< libdap::long_to_string(sizeFirst) << endl);
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - sizeSecond: "<< libdap::long_to_string(sizeSecond) << endl);

        if(sizeFirst < sizeSecond){
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - FNC Array first dimension is smaller than second." << endl);
            fncNodesDim = firstDim;
            fncFacesDim = secondDim;
        }
        else {
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - FNC Array second dimension is smaller than first." << endl);
            fncNodesDim = secondDim;
            fncFacesDim = firstDim;
        }

        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - fncFacesDim name: '" << fncArray->dimension_name(fncFacesDim) << "'" << endl);
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - fncFacesDim size: '" << fncArray->dimension_size(fncFacesDim,true) << "'" << endl);

        faceDimensionName = fncArray->dimension_name(fncFacesDim);
    }
    else {
        // There is already a faceDimensionName defined - possibly from loading face coordinate variables.

        // Does it match the name of the first or second dimensions of the fncArray? It better!
        if(faceDimensionName.compare(fncArray->dimension_name(firstDim)) == 0){
            fncNodesDim = secondDim;
            fncFacesDim = firstDim;
        }
        else if(faceDimensionName.compare(fncArray->dimension_name(secondDim)) == 0){
            fncNodesDim = firstDim;
            fncFacesDim = secondDim;
        }
        else {
            string msg = "The face coordinate dimension of the Face Node Connectivity variable '"+face_node_connectivity_var_name
                    +"' Has dimension name.'"+ fncArray->dimension_name(fncFacesDim) + "' which does not match the existing face coordinate dimension name '"
                    + faceDimensionName + "'";
            BESDEBUG("ugrid", msg << endl);
            throw Error(msg);
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Face dimension names match." << endl);

    }

    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Face dimension name: '" << faceDimensionName << "'" << endl);

    // Check to see if faceCount is initialized and do so if needed
    if(faceCount==0) {
        faceCount = fncArray->dimension_size(fncFacesDim,true);
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Face count: "<< libdap::long_to_string(faceCount) << endl);
    }
    else {
        // Make sure the face counts match.

        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Face count: "<< libdap::long_to_string(faceCount) << endl);
        if(faceCount!=fncArray->dimension_size(fncFacesDim,true) ){
            string msg = "The faces dimension of the Face Node Connectivity variable '"+face_node_connectivity_var_name
                    +"' Has size "+ libdap::long_to_string(fncArray->dimension_size(fncFacesDim,true))+ " which does not match the existing face count of "
                    + libdap::long_to_string(faceCount);
            BESDEBUG("ugrid", msg << endl);
            throw Error(msg);
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceNodeConnectivityArray() - Face counts match!" << endl);
    }

    faceNodeConnectivityArray = fncArray;

    BESDEBUG("ugrid", "TwoDMeshTopology::getFaceNodeConnectivityArray() - Got FCNA '"+fncArray->name()+"'" << endl);


}
/**
 * Returns the coordinate variables identified in the meshTopology variable's node_coordinates attribute.
 * throws an error if the node_coordinates attribute is missing, if the coordinates are not arrays, and
 * if the arrays are not all the same shape.
 */
void TwoDMeshTopology::ingestFaceCoordinateArrays(libdap::BaseType *meshTopology, libdap::DDS &dds)
{
    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - BEGIN Gathering face coordinate arrays..." << endl);

    string face_coordinates;
    AttrTable at = meshTopology->get_attr_table();

    AttrTable::Attr_iter iter_nodeCoors = at.simple_find(UGRID_FACE_COORDINATES);
    if (iter_nodeCoors != at.attr_end()) {
        face_coordinates = at.get_attr(iter_nodeCoors, 0);

        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Located '"<< UGRID_FACE_COORDINATES << "' attribute." << endl);


        if(faceCoordinateArrays==0)
            faceCoordinateArrays = new vector<libdap::Array *>();

        faceCoordinateArrays->clear();

        // Split the node_coordinates string up on spaces
        vector<string> faceCoordinateNames = split(face_coordinates, ' ');

        // Find each variable in the resulting list
        vector<string>::iterator coorName_it;
        for (coorName_it = faceCoordinateNames.begin();
                coorName_it != faceCoordinateNames.end(); ++coorName_it) {
            string faceCoordinateName = *coorName_it;

            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Processing face coordinate '"<< faceCoordinateName << "'." << endl);

            //Now that we have the name of the coordinate variable get it from the DDS!!
            BaseType *btp = dds.var(faceCoordinateName);
            if (btp == 0)
                throw Error(
                        "Could not locate the "  UGRID_FACE_COORDINATES  " variable named '"
                                + faceCoordinateName + "'! "
                                + "The mesh_topology variable is named "
                                + meshTopology->name());

            libdap::Array *newFaceCoordArray = dynamic_cast<libdap::Array*>(btp);
            if (newFaceCoordArray == 0) {
                throw Error(malformed_expr,
                        "Face coordinate variable '" + faceCoordinateName
                                + "' is not an Array type. It's an instance of "
                                + btp->type_name());
            }


            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Found face coordinate array '"<< faceCoordinateName << "'." << endl);

            // Coordinate arrays MUST be single dimensioned.
            if(newFaceCoordArray->dimensions(true) != 1){
                throw Error(malformed_expr,
                        "Face coordinate variable '" + faceCoordinateName
                                + "' has more than one dimension. That's just not allowed. It has "
                                + long_to_string(newFaceCoordArray->dimensions(true)) + " dimensions."
                               );
            }
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Face coordinate array '"<< faceCoordinateName << "' has a single dimension." << endl);


            // Make sure this node coordinate variable has the same size and name as all the others on the list - error if not true.
            string dimName = newFaceCoordArray->dimension_name(newFaceCoordArray->dim_begin());
            int    dimSize = newFaceCoordArray->dimension_size(newFaceCoordArray->dim_begin(), true);

            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - dimName: '"<< dimName << "' dimSize: "
                    << libdap::long_to_string(dimSize) << endl);

            if(faceDimensionName.empty()){
                faceDimensionName = dimName;
            }
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - faceDimensionName: '"<< faceDimensionName << "' " << endl);

            if(faceDimensionName.compare(dimName)!=0 ){
                throw Error(
                        "The face coordinate array '" + faceCoordinateName
                                + "' has the named dimension '"+ dimName + "' which differs from the expected  dimension name '"+faceDimensionName
                                +"'. The mesh_topology variable is named " + meshTopology->name());
            }
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Face dimension names match." << endl);


            if(faceCount == 0){
                faceCount = dimSize;
            }
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - faceCount: "<< libdap::long_to_string(faceCount) << endl);

            if(faceCount!=dimSize){
                throw Error(
                        "The face coordinate array '" + faceCoordinateName
                                + "' has a dimension size of " + libdap::long_to_string(dimSize) + " which differs from the the expected size of "
                                + libdap::long_to_string(faceCount) + " The mesh_topology variable is named "
                                + meshTopology->name());
            }
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Face counts match." << endl);

            // Add variable to faceCoordinateArrays.
            faceCoordinateArrays->push_back(newFaceCoordArray);
            BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Face coordinate array '"<< faceCoordinateName << "' ingested." << endl);
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - Located "<< libdap::long_to_string(faceCoordinateArrays->size()) << " face coordinate arrays." << endl);

    } else {
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - No Face Coordinates Found." << endl);
    }


    BESDEBUG("ugrid", "TwoDMeshTopology::ingestFaceCoordinateArrays() - DONE" << endl);

}


/**
 * Returns the coordinate variables identified in the meshTopology variable's node_coordinates attribute.
 * throws an error if the node_coordinates attribute is missing, if the coordinates are not arrays, and
 * if the arrays are not all the same shape.
 */
void TwoDMeshTopology::ingestNodeCoordinateArrays(libdap::BaseType *meshTopology, libdap::DDS &dds)
{
	BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - BEGIN Gathering node coordinate arrays..." << endl);

	string node_coordinates;
	AttrTable at = meshTopology->get_attr_table();

	AttrTable::Attr_iter iter_nodeCoors = at.simple_find(UGRID_NODE_COORDINATES);
	if (iter_nodeCoors != at.attr_end()) {
		node_coordinates = at.get_attr(iter_nodeCoors, 0);
	} else {
		throw Error(
				"Could not locate the " UGRID_NODE_COORDINATES " attribute in the " UGRID_MESH_TOPOLOGY
				" variable! The mesh_topology variable is named " + meshTopology->name());
	}
    BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Located '"<< UGRID_NODE_COORDINATES << "' attribute." << endl);

	if(nodeCoordinateArrays==0)
	    nodeCoordinateArrays = new vector<libdap::Array *>();

	nodeCoordinateArrays->clear();

	// Split the node_coordinates string up on spaces
	// TODO make this work on situations where multiple spaces in the node_coorindates string doesn't hose the split()
	vector<string> nodeCoordinateNames = split(node_coordinates, ' ');

	// Find each variable in the resulting list
	vector<string>::iterator coorName_it;
	for (coorName_it = nodeCoordinateNames.begin();
			coorName_it != nodeCoordinateNames.end(); ++coorName_it) {
		string nodeCoordinateName = *coorName_it;

	    BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Processing node coordinate '"<< nodeCoordinateName << "'." << endl);

		//Now that we have the name of the coordinate variable get it from the DDS!!
		BaseType *btp = dds.var(nodeCoordinateName);
		if (btp == 0)
			throw Error(
					"Could not locate the "  UGRID_NODE_COORDINATES  " variable named '"
					        + nodeCoordinateName + "'! "
							+ "The mesh_topology variable is named "
							+ meshTopology->name());

		libdap::Array *newNodeCoordArray = dynamic_cast<libdap::Array*>(btp);
		if (newNodeCoordArray == 0) {
			throw Error(malformed_expr,
					"Node coordinate variable '" + nodeCoordinateName
							+ "' is not an Array type. It's an instance of "
							+ btp->type_name());
		}
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Found node coordinate array '"<< nodeCoordinateName << "'." << endl);

		// Coordinate arrays MUST be single dimensioned.
		if(newNodeCoordArray->dimensions(true) != 1){
            throw Error(malformed_expr,
                    "Node coordinate variable '" + nodeCoordinateName
                            + "' has more than one dimension. That's just not allowed. It has "
                            + long_to_string(newNodeCoordArray->dimensions(true)) + " dimensions."
                           );
		}
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Node coordinate array '"<< nodeCoordinateName << "' has a single dimension." << endl);


		// Make sure this node coordinate variable has the same size and name as all the others on the list - error if not true.
        string dimName = newNodeCoordArray->dimension_name(newNodeCoordArray->dim_begin());
        int    dimSize = newNodeCoordArray->dimension_size(newNodeCoordArray->dim_begin(), true);

        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - dimName: '"<< dimName << "' dimSize: "
                << libdap::long_to_string(dimSize) << endl);

        if(nodeDimensionName.empty()){
            nodeDimensionName = dimName;
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - nodeDimensionName: '"<< nodeDimensionName << "' " << endl);
        if(nodeDimensionName.compare(dimName)!=0 ){
            throw Error(
                    "The node coordinate array '" + nodeCoordinateName
                            + "' has the named dimension '"+ dimName + "' which differs from the expected  dimension name '"+nodeDimensionName
                            +"'. The mesh_topology variable is named " + meshTopology->name());
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Node dimension names match." << endl);


        if(nodeCount == 0){
            nodeCount = dimSize;
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - nodeCount: "<< libdap::long_to_string(nodeCount) << endl);

        if(nodeCount!=dimSize){
            throw Error(
                    "The node coordinate array '" + nodeCoordinateName
                            + "' has a dimension size of " + libdap::long_to_string(dimSize) + " which differs from the the expected size of "
                            + libdap::long_to_string(nodeCount) + " The mesh_topology variable is named "
                            + meshTopology->name());
        }
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Node counts match." << endl);


		// Add variable to nodeCoordinateArrays vector.
		nodeCoordinateArrays->push_back(newNodeCoordArray);
        BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - Coordinate array '"<< nodeCoordinateName << "' ingested." << endl);


	}

    BESDEBUG("ugrid", "TwoDMeshTopology::ingestNodeCoordinateArrays() - DONE" << endl);

}



void TwoDMeshTopology::buildBasicGfTopology(){

    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Building GridFields objects for mesh_topology variable "<< myVar->name() << endl);
    // Start building the Grid for the GridField operation.
    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Constructing new GF::Grid for "<< name() << endl);
    gridTopology = new GF::Grid(name());

    // 1) Make the implicit nodes - same size as the node coordinate arrays
    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Building and adding implicit range Nodes to the GF::Grid" << endl);
    GF::AbstractCellArray *nodes = new GF::Implicit0Cells(nodeCount);
    // Attach the implicit nodes to the grid at rank 0
    gridTopology->setKCells(nodes, node);

    // @TODO Do I need to add implicit k-cells for faces (rank 2) if I plan to add range data on faces later?
    // Apparently not...


    // Attach the Mesh to the grid.
    // Get the face node connectivity cells (i think these correspond to the GridFields K cells of Rank 2)
    // FIXME Read this array once! It is read again below..
    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Building face node connectivity Cell array from the DAP version" << endl);
    GF::CellArray *faceNodeConnectivityCells = getFaceNodeConnectivityCells();


    // Attach the Mesh to the grid at rank 2
    // This 2 stands for rank 2, or faces.
    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Attaching Cell array to GF::Grid" << endl);
    gridTopology->setKCells(faceNodeConnectivityCells, face);

    // The Grid is complete. Now we make a GridField from the Grid
    BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Construct new GF::GridField from GF::Grid" << endl);
    d_inputGridField = new GF::GridField(gridTopology);
    // TODO Question for Bill: Can we delete the GF::Grid (tdmt->gridTopology) here?


    // We read and add the coordinate data (using GridField->addAttribute() to the GridField at
    // grid dimension/rank/dimension 0 (a.k.a. node)
    vector<libdap::Array *>::iterator ncit;
    for (ncit = nodeCoordinateArrays->begin(); ncit != nodeCoordinateArrays->end(); ++ncit) {
        libdap::Array *nca = *ncit;
        BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Adding node coordinate "<< nca->name() << " to GF::GridField at rank 0" << endl);
        GF::Array *gfa = extractGridFieldArray(nca,sharedIntArrays,sharedFloatArrays);
        d_inputGridField->AddAttribute(node, gfa);
    }

}



/**
 * @brief Builds a basic GF topology and then applies the passed filter restriction at the passed rank.
 */
void TwoDMeshTopology::buildRestrictedGfTopology(locationType loc, string filterExpression){
    buildBasicGfTopology();
    applyRestrictOperator(loc, filterExpression);
}


/**
 * @brief Builds the GridField Topology object and loads all all of the requested variables - OLDWAY
] */
void TwoDMeshTopology::buildGridFieldsTopology()
{

    buildBasicGfTopology();

	// For each range data variable associated with this MeshTopology read and add the  data to the GridField
	// At the appropriate rank.
	for (vector<MeshDataVariable *>::iterator mdv_it = rangeDataArrays->begin(); mdv_it != rangeDataArrays->end(); ++mdv_it) {
		MeshDataVariable *mdVar = *mdv_it;
		GF::Array *gfa = extractGridFieldArray(mdVar->getDapArray(),sharedIntArrays,sharedFloatArrays);
		BESDEBUG("ugrid", "TwoDMeshTopology::buildGridFieldsTopology() - Adding mesh data variable '"<< mdVar->getName() <<"' to GF::GridField at rank "<< mdVar->getGridLocation() << endl);
		d_inputGridField->AddAttribute(mdVar->getGridLocation(), gfa);
	}

}



/**
 * Takes a row major 3xN Face node connectivity DAP array
 * and converts it to a collection GF::Nodes organized as
 * 0,N,2N; 1,1+N,1+2N;
 *
 * This is the inverse operation to getGridFieldCellArrayAsDapArray()
 *
 *FIXME Make this use less memory. Certainly consider reading the values directly from
 *FIXME the DAP array (after it's read method has been called)
 */
GF::Node *TwoDMeshTopology::getFncArrayAsGFCells(libdap::Array *fncVar)
{

	BESDEBUG("ugrid", "TwoDMeshTopology::getFncArrayAsGFNodes() - BEGIN" << endl);


    int nodesPerFace = fncVar->dimension_size(fncNodesDim,true);
    int faceCount = fncVar->dimension_size(fncFacesDim,true);


	// interpret the array data as triangles
	GF::Node *cells = new GF::Node[ faceCount * nodesPerFace ];

	BESDEBUG("ugrid", "TwoDMeshTopology::getFncArrayAsGFNodes() - Reading DAP data into GF::Node array." << endl);
	GF::Node *temp_nodes = ugrid::extractArray<GF::Node>(fncVar);

	// Reorganize the cell ids so that cellids contains
	// the cells in three consecutive values (0,1,2; 3,4,5; ...).
	// The the values from  fncVar now in cellids2 and ar organized
	// as 0,N,2N; 1,1+N,1+2N; ...

	BESDEBUG("ugrid", "TwoDMeshTopology::getFncArrayAsGFNodes() - Re-packing and copying GF::Node array to result." << endl);
	for (int fIndex = 0; fIndex < faceCount; fIndex++) {
	    for(int nIndex=0; nIndex<nodesPerFace ; nIndex++){
	        cells[nodesPerFace * fIndex + nIndex] = temp_nodes[fIndex + (faceCount * nIndex)];

	    }
#if 0
		cellids[3 * fIndex] = cellids2[fIndex];
		cellids[3 * fIndex + 1] = cellids2[fIndex + faceCount];
		cellids[3 * fIndex + 2] = cellids2[fIndex + 2 * faceCount];
	}
#endif

    }

	BESDEBUG("ugrid", "TwoDMeshTopology::getFncArrayAsGFNodes() - Deleting intermediate GF::Node array." << endl);
	delete [] temp_nodes;

	BESDEBUG("ugrid", "TwoDMeshTopology::getFncArrayAsGFNodes() - DONE" << endl);

	return cells;
}

/**
 * Returns the value of the "start_index" attribute for the passed Array. If the start_index
 * is missing the value 0 is returned.
 */
int TwoDMeshTopology::getStartIndex(libdap::Array *array)
{
	AttrTable &at = array->get_attr_table();
	AttrTable::Attr_iter start_index_iter = at.simple_find(UGRID_START_INDEX);
	if (start_index_iter != at.attr_end()) {
		BESDEBUG("ugrid", "TwoDMeshTopology::getStartIndex() - Found the "<< UGRID_START_INDEX <<" attribute." << endl);
		AttrTable::entry *start_index_entry = *start_index_iter;
		if (start_index_entry->attr->size() == 1) {
			string val = (*start_index_entry->attr)[0];
			BESDEBUG("TwoDMeshTopology::getStartIndex", "value: " << val << endl);
			stringstream buffer(val);
			// what happens if string cannot be converted to an integer?
			int start_index;
			;

			buffer >> start_index;
			return start_index;
		} else {
			throw Error(malformed_expr,
					"Index origin attribute exists, but either no value supplied, or more than one value supplied.");
		}
	}
	return 0;
}

/**
 * Converts a row major 3xN Face node connectivity DAP array into a GF::CellArray
 */
GF::CellArray *TwoDMeshTopology::getFaceNodeConnectivityCells()
{

	BESDEBUG("ugrid", "TwoDMeshTopology::getFaceNodeConnectivityCells() - Building face node connectivity Cell " <<
			"array from the Array "<< faceNodeConnectivityArray->name() << endl);

	int nodesPerFace = faceNodeConnectivityArray->dimension_size(fncNodesDim);
	int total_size = nodesPerFace * faceCount;

	BESDEBUG("ugrid", "TwoDMeshTopology::getFaceNodeConnectivityCells() - Converting FNCArray to GF::Node array." << endl);
	sharedNodeArray = getFncArrayAsGFCells(faceNodeConnectivityArray);


	// adjust for the start_index (cardinal or ordinal array access)
	int startIndex = getStartIndex(faceNodeConnectivityArray);
	if (startIndex != 0) {
		BESDEBUG("ugrid", "TwoDMeshTopology::getFaceNodeConnectivityCells() - Applying startIndex to GF::Node array." << endl);
		for (int j = 0; j < total_size; j++) {
			sharedNodeArray[j] -= startIndex;
		}
	}
	// Create the cell array
	GF::CellArray *rankTwoCells = new GF::CellArray(sharedNodeArray, faceCount, nodesPerFace);


	BESDEBUG("ugrid", "TwoDMeshTopology::getFaceNodeConnectivityCells() - DONE" << endl);
	return rankTwoCells;

}



void TwoDMeshTopology::applyRestrictOperator(locationType loc, string filterExpression)
{
    // I think this function could be done with just the following single line:
    // resultGridField = GF::RefRestrictOp::Restrict(filterExpression,loc,d_inputGridField);

	// Build the restriction operator;
	BESDEBUG("ugrid", "TwoDMeshTopology::applyRestrictOperator() - Constructing new GF::RestrictOp using user "<<
			"supplied 'dimension' value and filter expression combined with the GF:GridField " << endl);
	GF::RestrictOp op = GF::RestrictOp(filterExpression, loc, d_inputGridField);


	// Apply the operator and get the result;
	BESDEBUG("ugrid", "TwoDMeshTopology::applyRestrictOperator() - Applying GridField operator." << endl);
	GF::GridField *resultGF = op.getResult();
	resultGridField = resultGF;
}


/**
 *
 */
static void rDAWorker(
        libdap::Array *dapArray,
        libdap::Array::Dim_iter thisDim,
        libdap::Array::Dim_iter locationCoordinateDim,
        locationType gridLocation,
        NDimensionalArray *results) {

    if(thisDim==dapArray->dim_end()){

       // GF::Array *gfa = extractGridFieldArray(dapArray,sharedIntArrays,sharedFloatArrays);

        //resultGridField->Bind(gridLocation,gfa);
        //resultGridField->GetGrid()->normalize();

        //libdap::Array *resultRangeVar = getGFAttributeAsDapArray( dapArray, gridLocation, resultGridField);

        //resultGridField->unBind(gridLocation,gfa);


    }
    else {
        libdap::Array::Dim_iter nextDim = thisDim;
        nextDim++; // now it's the next dimension.

        if(thisDim==locationCoordinateDim){

            if(nextDim!=dapArray->dim_end()){
                string msg = "rDAWorker() - The location coordinate dimension is not the last dimension in the array. Hyperslab subsetting of this dimension is not supported.";
                BESDEBUG("ugrid", msg << endl);
                throw Error(malformed_expr, msg);

            }


            rDAWorker(dapArray, nextDim, locationCoordinateDim, gridLocation, results);
        }
        else {


            unsigned int start  = dapArray->dimension_start(thisDim, true);
            unsigned int stride = dapArray->dimension_stride(thisDim, true);
            unsigned int stop   = dapArray->dimension_stop(thisDim, true);

            for(unsigned int dimIndex=start; dimIndex<=stop ; dimIndex+=stride){
                dapArray->add_constraint(thisDim,dimIndex,1,dimIndex);
                rDAWorker(dapArray, nextDim, locationCoordinateDim, gridLocation, results);
            }

            // Reset the constraint for this dimension.
            dapArray->add_constraint(thisDim,start,stride,stop);
        }

    }


}



/**
 * Now, for each variable array on the mesh we have to hyper-slab the array such that all the dimensions, with the
 * exception of the rank/location dimension (the one that is either the number of nodes, edges, or faces depending
 * on which of these locations is the data is associated with), have size one.
 *  Each of these slabs is then added to the GridField, subset and the result must be packed back
 *  into the result array so that things work out.
**/
libdap::Array *TwoDMeshTopology::restrictDapArray(libdap::Array *dapArray, libdap::Array::Dim_iter locationCoordinateDim, locationType gridLocation){


    // We want the manipulate the Array's Dimensions so that only a single dimensioned slab of the location coordinate dimension
    // is read at a time. We need to cache the original constrained dimensions so that we can build the correct collection of
    // location coordinate dimension slabs.

    // What's the shape of the requested range variable array?
    vector<unsigned int> arrayShape(dapArray->dimensions(true));
    NDimensionalArray::computeConstrainedShape(dapArray, &arrayShape );

    // Now,
    arrayShape[dapArray->dimensions(true) - 1] = resultGridField->GetGrid()->countKCells(node);

    NDimensionalArray *result = new NDimensionalArray(&arrayShape, dapArray->var()->type());

    rDAWorker(dapArray, dapArray->dim_begin(), locationCoordinateDim, gridLocation, result);

    libdap::Array  *myResult = 0;

    return myResult;

}



void TwoDMeshTopology::restrictRange( vector<BaseType *> *results){


     BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - BEGIN" << endl);


     BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - Normalizing Grid." << endl);
     resultGridField->GetGrid()->normalize();

     // Add the coordinate node arrays to the response.
     BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - Adding the coordinate node arrays to the response." << endl);
     vector<libdap::Array *>::iterator it;
     for (it = nodeCoordinateArrays->begin(); it != nodeCoordinateArrays->end(); ++it) {
         libdap::Array *sourceCoordinateArray = *it;
         libdap::Array *resultCoordinateArray = getGFAttributeAsDapArray(sourceCoordinateArray, node, resultGridField );
         results->push_back(resultCoordinateArray);
     }
    // Add the new face node connectivity array - make sure it has the same attributes as the original.
    BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - Adding the new face node connectivity array to the response." << endl);
    libdap::Array *resultFaceNodeConnectivityDapArray = getGridFieldCellArrayAsDapArray(resultGridField, faceNodeConnectivityArray);
    results->push_back(resultFaceNodeConnectivityDapArray);


    // Add the range variable data arrays to the response.
    BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - Processing range variable data arrays." << endl);
    vector<MeshDataVariable *>::iterator mdvIt;
    for (mdvIt = rangeDataArrays->begin(); mdvIt != rangeDataArrays->end(); ++mdvIt) {
        MeshDataVariable *mdv = *mdvIt;


        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Processing MeshDataVariable '"<< mdv->getName() << "' array type is "<< mdv->getDapArray()->type_name() << endl);

        libdap::Array *resultRangeVar = restrictDapArray(mdv->getDapArray(), mdv->getLocationCoordinateDimension(), mdv->getGridLocation());

        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Retrieved restricted DAP Array for MeshDataVariable '"<< mdv->getName() << "'" << endl);
        results->push_back(resultRangeVar);
        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - DAP Array  added to result" << endl);

        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - DAP Array  added to result" << endl);
    }






    BESDEBUG("ugrid", "TwoDMeshTopology::restrictRange() - END" << endl);

}








/**
 * Builds the DAP response content from the GF::GridField result object.
 */
void TwoDMeshTopology::convertResultGridFieldToDapObjects(vector<BaseType *> *results)
{

	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - BEGIN" << endl);

	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Normalizing Grid." << endl);
	resultGridField->GetGrid()->normalize();

	// Add the coordinate node arrays to the response.
	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Adding the coordinate node arrays to the response." << endl);
	vector<libdap::Array *>::iterator it;
	for (it = nodeCoordinateArrays->begin(); it != nodeCoordinateArrays->end(); ++it) {
		libdap::Array *sourceCoordinateArray = *it;
		libdap::Array *resultCoordinateArray = getGFAttributeAsDapArray(sourceCoordinateArray, node, resultGridField );
		results->push_back(resultCoordinateArray);
	}

	// Add the range variable data arrays to the response.
	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Adding the range variable data arrays to the response." << endl);
	vector<MeshDataVariable *>::iterator mdvIt;
	for (mdvIt = rangeDataArrays->begin(); mdvIt != rangeDataArrays->end(); ++mdvIt) {
		MeshDataVariable *mdv = *mdvIt;
	    BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Processing MeshDataVariable '"<< mdv->getName() << "' array type is "<< mdv->getDapArray()->type_name() << endl);
		libdap::Array *resultRangeVar = getGFAttributeAsDapArray( mdv->getDapArray(), mdv->getGridLocation(), resultGridField);
        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Retrieved restricted DAP Array for MeshDataVariable '"<< mdv->getName() << "'" << endl);
		results->push_back(resultRangeVar);
        BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - DAP Array  added to result" << endl);
	}

	// Add the new face node connectivity array - make sure it has the same attributes as the original.
	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - Adding the new face node connectivity array to the response." << endl);
	libdap::Array *resultFaceNodeConnectivityDapArray = getGridFieldCellArrayAsDapArray(resultGridField, faceNodeConnectivityArray);
	results->push_back(resultFaceNodeConnectivityDapArray);

	BESDEBUG("ugrid", "TwoDMeshTopology::convertResultGridFieldToDapObjects() - END" << endl);

}




/**
 * Get a new 3xN DAP Array of Int32 with the same name, attributes, and dimension names
 * as the templateArray. Make the new array's second dimension size N.
 * Returns a DAP Array with an Int32 type template.
 */
libdap::Array *TwoDMeshTopology::getNewFncDapArray(libdap::Array *templateArray, int N)
{

	// Is the template array a 2D array?
	int dimCount = templateArray->dimensions(true);
	if (dimCount != 2)
		throw Error(
				"Expected a 2 dimensional array. The array '"
						+ templateArray->name() + "' has "
						+ libdap::long_to_string(dimCount) + " dimensions.");

	// Is the template array really 3xN?
	libdap::Array::Dim_iter di = templateArray->dim_begin();
	if (di->c_size != 3) {
		string msg =
				"Expected a 2 dimensional array with shape of 3xN! The array "
						+ templateArray->name() + " has a first "
						+ "dimension of size " + libdap::long_to_string(di->c_size);
		BESDEBUG("ugrid", "TwoDMeshTopology::getNewFcnDapArray() - "<< msg << endl);
		throw Error(malformed_expr, msg);
	}

	// Get a new template variable for our new array (should be just like the template for the source array)
	//BaseType *arrayTemplate = getDapVariableInstance(templateArray->var(0)->name(),templateArray->var(0)->type());
	libdap::Array *newArray = new libdap::Array(templateArray->name(),
			new Int32(templateArray->name()));

	//Add the first dimension (size 3 same same as template array's first dimension)
	newArray->append_dim(3, di->name);

	// Add the second dimension to the result array, but use only the name from the template array's
	// second dimension. The size will be from the passed parameter N
	di++;
	newArray->append_dim(N, di->name);

	newArray->set_attr_table(templateArray->get_attr_table());

	// make the new array big enough to hold all the values.
	newArray->reserve_value_capacity(3 * N);

#if 0
	cerr<<"getNewFcnDapArray() -"<<endl<<endl;
	cerr << "Newly minted Array: "<< endl;
	newArray->print_val(cerr);
	cerr<<endl<<endl;

#endif

	return newArray;

}

/**
 * Takes a GF::GridField, extracts it's rank2 GF::CellArray. The GF::CellArray content is
 * extracted and re-packed into a 3xN DAP Array. This is the inverse operation to
 * getFncArrayAsGFNodes()
 */
libdap::Array *TwoDMeshTopology::getGridFieldCellArrayAsDapArray(GF::GridField *resultGridField, libdap::Array *sourceFcnArray)
{

	BESDEBUG("ugrid", "TwoDMeshTopology::getGridFieldCellArrayAsDapArray() - BEGIN" << endl);

	// Get the rank 2 k-cells from the GridField object.
	GF::CellArray* gfCellArray = (GF::CellArray*) (resultGridField->GetGrid()->getKCells(2));

	// This is a vector of size N holding vectors of size 3
	vector<vector<int> > nodes2 = gfCellArray->makeArrayInts();

	libdap::Array *resultFncDapArray = getNewFncDapArray(sourceFcnArray, nodes2.size());

	// Make a vector to hold the re-packed cell nodes.
	vector<dods_int32> rowMajorNodes;

	// Re-pack the mesh nodes.
	for (unsigned int firstDim = 0; firstDim < 3; firstDim++) {
		for (unsigned int secondDim = 0; secondDim < nodes2.size();
				secondDim++) {
			dods_int32 val = nodes2.at(secondDim).at(firstDim);
			rowMajorNodes.push_back(val);
		}
	}

#if 0    // This is debugging - keep it around we might need it.
	cerr << "getGridFieldCellArrayAsDapArray() - rowMajorNodes: " << endl << "{";
	for (unsigned int j=0; j < rowMajorNodes.size(); j++) {
		dods_int32 val = rowMajorNodes.at(j);
		cerr << val << ", ";
	}
	cerr << "}" << endl;
#endif

	// Add them to the DAP array.
	resultFncDapArray->set_value(rowMajorNodes, rowMajorNodes.size());

#if 0    // This is debugging - keep it around we might need it.
	cerr << "getGridFieldCellArrayAsDapArray() - DAP Array: "<< endl;
	resultFcnDapArray->print_val(cerr);
#endif
	BESDEBUG("ugrid", "TwoDMeshTopology::getGridFieldCellArrayAsDapArray() - DONE" << endl);

	return resultFncDapArray;

}



static libdap::Array::Dim_iter copySizeOneDimensions(libdap::Array *sourceArray, libdap::Array *dapArray){

    libdap::Array::Dim_iter dataDimension;
    for (libdap::Array::Dim_iter srcArrIt = sourceArray->dim_begin();   srcArrIt != sourceArray->dim_end() ;++srcArrIt) {

       // Get the original dimension size
        int dimSize = sourceArray->dimension_size(srcArrIt,true);
        string dimName = sourceArray->dimension_name(srcArrIt);

        // Preserve single dimensions
        if(dimSize==1){
            BESDEBUG("ugrid", "TwoDMeshTopology::copySizeOneDimensions() - Adding size one dimension '"<< dimName << "' from source array into result." << endl);
            dapArray->append_dim(dimSize, dimName);
        }
        else {
            BESDEBUG("ugrid", "TwoDMeshTopology::copySizeOneDimensions() - Located data dimension '"<< dimName << "' in source array." << endl);

            dataDimension = srcArrIt;
        }
    }
    BESDEBUG("ugrid", "TwoDMeshTopology::copySizeOneDimensions() - Returning dimension iterator pointing to '"<< sourceArray->dimension_name(dataDimension) << "'." << endl);
    return dataDimension;
}



/**
 * Retrieves a single dimensional rank 0 GF attribute array from a GF::GridField and places the data into
 * DAP array of the appropriate type.
 */
libdap::Array *TwoDMeshTopology::getGFAttributeAsDapArray(libdap::Array *sourceArray, locationType rank, GF::GridField *resultGridField)
{

	BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - BEGIN" << endl);

	// The result variable is assumed to be bound to the GridField at 'rank'
	// Try to get the Attribute from 'rank' with the same name as the source array
	BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - Retrieving GF::GridField Attribute '" <<
			sourceArray->name() << "'" << endl);
	GF::Array* gfa = resultGridField->GetAttribute(rank, sourceArray->name());

	libdap::Array *dapArray;
	BaseType *templateVar = sourceArray->var();
	string dimName;

	switch (templateVar->type()) {
	case dods_byte_c:
	case dods_uint16_c:
	case dods_int16_c:
	case dods_uint32_c:
	case dods_int32_c: {
		// Get the data
		BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - GF::Array was made from some type of int, retrieve it as such." << endl);
		vector<dods_int32> GF_ints = gfa->makeArray();
		// Make a DAP array to put the data into.
		dapArray = new libdap::Array(sourceArray->name(), new Int32(sourceArray->name()));

		// copy the dimensions whose size is "1" from the source array to the result.
		libdap::Array::Dim_iter dataDimension = copySizeOneDimensions(sourceArray, dapArray);

		// Add the result dimension
		dimName = sourceArray->dimension_name(dataDimension);
	    BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - Adding data dimension: '"<< dimName << "'" << endl);

		dapArray->append_dim(GF_ints.size(), dimName);

		// Add the data
		dapArray->set_value(GF_ints, GF_ints.size());
		break;
	}
	case dods_float32_c:
	case dods_float64_c: {
		// Get the data
		BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - GF::Array was made from some type of float, retrieve it as such." << endl);
		vector<dods_float64> GF_floats = gfa->makeArrayf();
		// Make a DAP array to put the data into.
		dapArray = new libdap::Array(sourceArray->name(), new Float64(sourceArray->name()));

        // copy the dimensions whose size is "1" from the source array to the result.
        libdap::Array::Dim_iter dataDimension = copySizeOneDimensions(sourceArray, dapArray);

        // Add the result dimension
        dimName = sourceArray->dimension_name(dataDimension);
        dapArray->append_dim(GF_floats.size(), dimName);

		// Add the data
		dapArray->set_value(GF_floats, GF_floats.size());
		break;
	}
	default:
		throw InternalErr(__FILE__, __LINE__,
				"Unknown DAP type encountered when converting to gridfields array");
	}

	// Copy the source objects attributes.
	dapArray->set_attr_table(sourceArray->get_attr_table());

	BESDEBUG("ugrid", "TwoDMeshTopology::getGFAttributeAsDapArray() - DONE" << endl);

	return dapArray;
}








} // namespace ugrid
