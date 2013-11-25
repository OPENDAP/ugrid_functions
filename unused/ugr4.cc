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

// NOTE: This file is built only when the gridfields library is linked with
// the netcdf_handler (i.e., the handler's build is configured using the
// --with-gridfields=... option to the 'configure' script).

#include "config.h"

#include <limits.h>

#include <cstdlib>      // used by strtod()
#include <cerrno>
#include <cmath>
#include <iostream>
#include <sstream>
//#include <cxxabi.h>

#define DODS_DEBUG

#include "BaseType.h"
#include "Int32.h"
#include "Str.h"
#include "Array.h"
#include "Structure.h"
#include "Error.h"

#include "BESDebug.h"
#include "util.h"

#include "ugrid_utils.h"
#include "MeshDataVariable.h"
#include "TwoDMeshTopology.h"
#include "NDimensionalArray.h"

#include "ugr4.h"


using namespace std;
using namespace libdap;

namespace ugrid {

/**
 * Function syntax
 */
static string ugrSyntax =
		"ugr4(dim:int32, rangeVariable:string, [rangeVariable:string, ... ] condition:string)";

/**
 * Function Arguments
 */
struct UgridRestrictArgs {
	locationType dimension;
	vector<libdap::Array *> rangeVars;
	string filterExpression;
};


/**
 * Evaluates the rangeVar and determines which meshTopology it is associated with. If one hasn't been found
 * a new mesh topology is created. Once the associated mesh topology had been found (or created), the rangeVar
 * is added to the vector of rangeVars held by the mesh topology for later evaluation.
 */
static void addRangeVar(DDS *dds, libdap::Array *rangeVar, map<string, vector<MeshDataVariable *> *> *rangeVariables) {

	MeshDataVariable *mdv = new MeshDataVariable();
	mdv->init(rangeVar);
	string meshVarName = mdv->getMeshName();

    BaseType *meshVar = dds->var(meshVarName);

    if(meshVar == 0){
        string msg = "The range variable '"+mdv->getName()+"' references the mesh variable '"+meshVarName+
                "' which cannot be located in this dataset.";
        BESDEBUG("ugrid", "addRangeVar() - " << msg  << endl);
        throw new Error(no_such_variable,msg);
    }


	// Get the rangeVariable vector for this mesh name from the map.
	vector<MeshDataVariable *> *requestedRangeVarsForMesh;
	map<string, vector<MeshDataVariable *> *>::iterator mit = rangeVariables->find(meshVarName);
	if(mit == rangeVariables->end()){
		// Not there? Make a new one.
		BESDEBUG("ugrid", "addRangeVar() - MeshTopology object for '" << meshVarName <<"' does NOT exist. Getting a 'new' one... "  << endl);

		requestedRangeVarsForMesh =  new vector<MeshDataVariable *>();
		(*rangeVariables)[meshVarName] = requestedRangeVarsForMesh;
	}
	else {
		// Sweet! Found it....
		BESDEBUG("ugrid", "addRangeVar() - MeshTopology object for '" << meshVarName <<"' exists. Retrieving... "  << endl);
		requestedRangeVarsForMesh = mit->second;
	}

	requestedRangeVarsForMesh->push_back(mdv);

}

/**
 * Process the functions arguments and return the structure containing their values.
 */
static UgridRestrictArgs processUgrArgs(int argc, BaseType *argv[]) {

	BESDEBUG( "ugrid", "processUgrArgs() - BEGIN" << endl);

	UgridRestrictArgs args;
	args.rangeVars = vector<libdap::Array *>();

	// Check number of arguments; BESDEBUG is a macro. Use #define
	// DODS_DEBUG to activate the debugging stuff.
	if (argc < 3)
		throw Error(malformed_expr,
				"Wrong number of arguments to ugrid restrict function: "
						+ ugrSyntax + " was passed "
						+ long_to_string(argc) + " argument(s)");

	BaseType * bt;

	// ---------------------------------------------
	// Process the first arg, which is the rank of the Restriction Clause
	bt = argv[0];
	if (bt->type() != dods_int32_c)
		throw Error(malformed_expr,
				"Wrong type for first argument, expected DAP Int32. "
						+ ugrSyntax + "  was passed a/an "
						+ bt->type_name());
	args.dimension = (locationType) dynamic_cast<Int32&>(*argv[0]).value();


	// ---------------------------------------------
	// Process the last argument, the relational expression used to restrict the ugrid content.
	bt = argv[argc - 1];
	if (bt->type() != dods_str_c)
		throw Error(malformed_expr,
				"Wrong type for third argument, expected DAP String. "
						+ ugrSyntax + "  was passed a/an "
						+ bt->type_name());
	args.filterExpression = dynamic_cast<Str&>(*bt).value();


	// --------------------------------------------------
	// Process the range variables selected by the user.
	// We know that argc>=3, because we checked so the
	// following loop will try to find at least one rangeVar,
	// and it won't try to process the first or last members
	// of argv.
	for (int i = 1; i < (argc - 1); i++) {
		bt = argv[i];
		if (bt->type() != dods_array_c)
			throw Error(malformed_expr,
					"Wrong type for second argument, expected DAP Array. "
							+ ugrSyntax + "  was passed a/an "
							+ bt->type_name());

		libdap::Array *newRangeVar = dynamic_cast<libdap::Array*>(bt);
		if (newRangeVar == 0) {
			throw Error(malformed_expr,
					"Wrong type for second argument. " + ugrSyntax
							+ "  was passed a/an " + bt->type_name());
		}
		args.rangeVars.push_back(newRangeVar);
	}

	BESDEBUG("ugrid", "processUgrArgs() - END" << endl);

	return args;

}


static string arrayState(libdap::Array *dapArray, string indent){

    libdap::Array::Dim_iter thisDim;
    stringstream arrayState;
    arrayState << endl;
    arrayState << indent << "dapArray: '" << dapArray->name() << "'   " ;
    arrayState << "type: '" << dapArray->type_name( ) << "'  " ;

    arrayState << "read(): " << dapArray->read() << "  ";
    arrayState << "read_p(): " << dapArray->read_p() << "  ";
    arrayState << endl;

    for(thisDim= dapArray->dim_begin(); thisDim!=dapArray->dim_end() ;thisDim++){

    arrayState << indent << "  thisDim:  '" << dapArray->dimension_name(thisDim) << "' ";
    arrayState << indent << "    start:  " << dapArray->dimension_start(thisDim)  << " ";
    arrayState << indent << "    stride: " << dapArray->dimension_stride(thisDim)  << " ";
    arrayState << indent << "    stop:   " << dapArray->dimension_stop(thisDim)  << " ";
    arrayState << endl;
    }
    return arrayState.str();
}



/**
 *
 */
static void rDAWorker(
        MeshDataVariable *mdv,
        libdap::Array::Dim_iter thisDim,
        DDS *dds,
        locationType restrictLocation,
        string filterExpression,
        NDimensionalArray *results) {

    libdap::Array *dapArray = mdv->getDapArray();

    if(thisDim==dapArray->dim_end()){



        BESDEBUG("ugrid", "rdaWorker() - Arrived At Last Dimension"<<endl);

        dapArray->set_read_p(false);

        BESDEBUG("ugrid", "rdaWorker() - dap array: "<< arrayState(dapArray, "    "));

        TwoDMeshTopology *tdmt = new TwoDMeshTopology();
        tdmt->init(mdv->getMeshName(), dds);
        tdmt->addDataVariable(mdv);
        tdmt->buildGridFieldsTopology();
        tdmt->applyRestrictOperator(restrictLocation, filterExpression);

        vector<unsigned int> lastDimHyperSlabLocation;
        NDimensionalArray::retrieveLastDimHyperSlabLocationFromConstrainedArrray(dapArray,&lastDimHyperSlabLocation);

        BESDEBUG("ugrid", "rdaWorker() - lastDimHyperSlabLocation: "<< NDimensionalArray::vectorToIndices(&lastDimHyperSlabLocation) << endl);

        unsigned int elementCount;

        void *slab;
        results->getLastDimensionHyperSlab(&lastDimHyperSlabLocation,&slab,&elementCount);


        tdmt->getResultGFAttributeValues(dapArray,mdv->getGridLocation(),slab);

        //tdmt->getResultGFAttributeValues(dapArray->name(),dapArray->var()->type(),mdv->getGridLocation(),slab);



        stringstream s;
        for(int i=0; i<elementCount ; i++){
            double f = ((double *) slab)[i];
            s << "slab[" << i << "]: " << f << endl;
        }
        BESDEBUG("ugrid", "rdaWorker() - Retrieved slab: "<< endl << s.str());

        BESDEBUG("ugrid", "rdaWorker() - results: " << results->toString());


        delete tdmt;

    }
    else {
        libdap::Array::Dim_iter locationCoordinateDim = mdv->getLocationCoordinateDimension();

        libdap::Array::Dim_iter nextDim = thisDim;
        nextDim++; // now it's the next dimension.
        BESDEBUG("ugrid", "rdaWorker() - thisDim: '" << dapArray->dimension_name(thisDim) << "'" << endl);
        BESDEBUG("ugrid", "rdaWorker() - locationCoordinateDim: '" << dapArray->dimension_name(locationCoordinateDim) << "'" << endl);

        if(thisDim==locationCoordinateDim){
            BESDEBUG("ugrid", "rdaWorker() - Found locationCoordinateDim" << endl);

            if(nextDim!=dapArray->dim_end()){
                string msg = "rDAWorker() - The location coordinate dimension is not the last dimension in the array. Hyperslab subsetting of this dimension is not supported.";
                BESDEBUG("ugrid", msg << endl);
                throw Error(malformed_expr, msg);

            }

            rDAWorker(mdv, nextDim, dds,restrictLocation, filterExpression, results);
        }
        else {


            unsigned int start  = dapArray->dimension_start(thisDim, true);
            unsigned int stride = dapArray->dimension_stride(thisDim, true);
            unsigned int stop   = dapArray->dimension_stop(thisDim, true);

            BESDEBUG("ugrid", "rdaWorker() - array state: " << arrayState(dapArray, "    "));

            for(unsigned int dimIndex=start; dimIndex<=stop ; dimIndex+=stride){
                dapArray->add_constraint(thisDim,dimIndex,1,dimIndex);
                rDAWorker(mdv, nextDim, dds,restrictLocation, filterExpression, results);
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
static libdap::Array *restrictRangeVariableByOneDHyperSlab(
        MeshDataVariable *mdv,
        long restrictedSlabSize,
        string meshVariableName,
        DDS *dds,
        locationType restrictLocation,
        string filterExpression
){

    libdap:Array *sourceDapArray = mdv->getDapArray();

    BESDEBUG("ugrid", "restrictDapArrayByOneDHyperSlab() - locationCoordinateDim: '" <<
            sourceDapArray->dimension_name(mdv->getLocationCoordinateDimension()) << "'" << endl);




    // We want the manipulate the Array's Dimensions so that only a single dimensioned slab of the location coordinate dimension
    // is read at a time. We need to cache the original constrained dimensions so that we can build the correct collection of
    // location coordinate dimension slabs.

    // What's the shape of the requested range variable array?
    vector<unsigned int> arrayShape(sourceDapArray->dimensions(true));
    NDimensionalArray::computeConstrainedShape(sourceDapArray, &arrayShape );

    stringstream msg;
    for(int i=0; i< arrayShape.size(); i++){
        msg << "[" << arrayShape[i] << "]";
    }
    BESDEBUG("ugrid", "restrictDapArrayByOneDHyperSlab() - arrayShape" << msg.str() << endl);
    msg.str(std::string());

    // Now, we know that the restricted slab size the new (restricted)n size of the last dimension,
    // so we make the shape reflect that

    arrayShape[sourceDapArray->dimensions(true) - 1] = restrictedSlabSize;
    libdap::Type dapType = sourceDapArray->var()->type();


    BESDEBUG("ugrid", "restrictDapArrayByOneDHyperSlab() - Restricted HyperSlab has  "<< restrictedSlabSize << " elements." << endl);
    BESDEBUG("ugrid", "restrictDapArrayByOneDHyperSlab() - Array is of type '"<< libdap::type_name(dapType) << "'" << endl);

    for(int i=0; i< arrayShape.size(); i++){
        msg << "[" << arrayShape[i] << "]";
    }
    BESDEBUG("ugrid", "restrictDapArrayByOneDHyperSlab() - arrayShape" << msg.str() << endl);
    msg.str(std::string());

    Type resultType = ugrid::getGridfieldsReturnType(dapType);

    NDimensionalArray *result = new NDimensionalArray(&arrayShape, resultType);

    //rDAWorker(dapArray, dapArray->dim_begin(), locationCoordinateDim, gridLocation, result);
    rDAWorker(mdv, sourceDapArray->dim_begin(), dds, restrictLocation, filterExpression, result);


    libdap::Array  *resultDapArray = result->getArray(sourceDapArray);

    delete result;

    return resultDapArray;

}






/**
 Subset an irregular mesh (aka unstructured grid).

 @param argc Count of the function's arguments
 @param argv Array of pointers to the functions arguments
 @param dds Reference to the DDS object for the complete dataset.
 This holds pointers to all of the variables and attributes in the
 dataset.
 @param btpp Return the function result in an instance of BaseType
 referenced by this pointer to a pointer. We could have used a
 BaseType reference, instead of pointer to a pointer, but we didn't.
 This is a value-result parameter.

 @return void

 @exception Error Thrown If the Array is not a one dimensional
 array. */
void ugr4(int argc, BaseType *argv[], DDS &dds, BaseType **btpp)
{
	BESDEBUG("ugrid", "ugr4() - BEGIN" << endl);

	string info = string("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
			+ "<function name=\"ugr4\" version=\"0.1\">\n"
			+ "Server function for Unstructured grid operations.\n" + "usage: "
			+ ugrSyntax + "\n"+ "</function>";

	if (argc == 0) {
		Str *response = new Str("info");
		response->set_value(info);
		*btpp = response;
		return;
	}

	// Process and QC the arguments
	UgridRestrictArgs args = processUgrArgs(argc, argv);

	// Each range variable is associated with a "mesh" i.e. a mesh topology variable. Since there may be more than one mesh in a
	// dataset, and the user may request more than one range variable for each mesh we need to sift through the list of requested
	// range variables and organize them by mesh topology variable name.
	map<string, vector<MeshDataVariable *> *> *meshToRangeVarsMap = new map<string, vector<MeshDataVariable *> *>();

	// For every Range variable in the arguments list, locate it and ingest it.
	vector<libdap::Array *>::iterator it;
	for (it = args.rangeVars.begin(); it != args.rangeVars.end(); ++it) {
	    addRangeVar(&dds, *it, meshToRangeVarsMap);
	}
	BESDEBUG("ugrid", "ugr4() - The user requested "<< args.rangeVars.size() << " range data variables." << endl);
	BESDEBUG("ugrid", "ugr4() - The user's request referenced "<< meshToRangeVarsMap->size() << " mesh topology variables." << endl);

	// ----------------------------------
	// OK, so up to this point we have not read any data from the data set, but we have QC'd the inputs and verified that
	// it looks like the request is consistent with the semantics of the dataset.
	// Now it's time to read some data and pack it into the GridFields library...


    // TODO This returns a single structure but it would make better sense to the
    // world if it could return a vector of objects and have them appear at the
    // top level of the DDS.
    // FIXME fix the names of the variables in the mesh_topology attributes
    // If the server side function can be made to return a DDS or a collection of BaseType's then the
    // names won't change and the original mesh_topology variable and it's metadata will be valid
	Structure *dapResult = new Structure("ugr_result");


	// Since we only want each ugrid structure to appear in the results one time  (cause otherwise we might be trying to add
	// the same variables with the same names to the result multiple times.) we grind on this by iterating over the
	// names of the mesh topology names.
    vector<MeshDataVariable *>::iterator rvit;
    map<string, vector<MeshDataVariable *> *>::iterator mit;
	for (mit = meshToRangeVarsMap->begin(); mit != meshToRangeVarsMap->end(); ++mit) {

		string meshVariableName = mit->first;
		vector<MeshDataVariable *> *requestedRangeVarsForMesh = mit->second;

		// When we built the meshToRangeVarsMap we QC'd this so we already know that
		// the meshVar exists - no need to re-check.
        BaseType *meshVar = dds.var(meshVariableName);

        vector<BaseType *> dapResults;

        // Building the restricted TwoDMeshTopology without adding any range variables and then converting the result
        // Grid field to Dap Objects should return all of the Ugrid structural stuff - mesh variable, node coordinate variables,
        // face and edge coordinate variables if present.
        BESDEBUG("ugrid", "ugr4() - Adding restricted mesh_topology structure for mesh '" << meshVariableName << "' to DAP response." << endl);

        TwoDMeshTopology *tdmt = new TwoDMeshTopology();
        tdmt->init(meshVariableName, &dds);
        tdmt->buildRestrictedGfTopology(args.dimension, args.filterExpression);
        tdmt->convertResultGridFieldStructureToDapObjects(&dapResults);


        BESDEBUG("ugrid", "ugr4() - Restriction of mesh_topology '"<< tdmt->getMeshVariable()->name() << "' structure completed." << endl);

        // now that we have the mesh topology variable we are going to look at each of the requested
        // range variables (aka MeshDataVariable instances) and we're going to subset that using the
        // gridfields library and add its subset version to the results.
	    for(rvit=requestedRangeVarsForMesh->begin(); rvit!=requestedRangeVarsForMesh->end(); rvit++){
	        MeshDataVariable *mdv = *rvit;

	        BESDEBUG("ugrid", "ugr4() - Processing MeshDataVariable  '"<< mdv->getName() << "' ." << endl);

	        tdmt->setLocationCoordinateDimension(mdv);

	        /**
	         * Here is where we will do the range variable sub-setting including decomposing the requested variable
	         * into 1-dimensional hyper-slabs that can be fed into the gridfields library
	         */
            libdap:Array *restrictedRangeVarArray =
                    restrictRangeVariableByOneDHyperSlab(mdv, tdmt->getResultGridSize(mdv->getGridLocation()), meshVariableName, &dds, args.dimension, args.filterExpression);

            BESDEBUG("ugrid", "ugr4() - Adding resulting dapArray  '"<< restrictedRangeVarArray->name() << "' to dapResults." << endl);

            dapResults.push_back(restrictedRangeVarArray);


	    }

        delete tdmt;


		BESDEBUG("ugrid", "ugr4() - Adding GF::GridField results to DAP structure " << dapResult->name() << endl);
		for (vector<BaseType *>::iterator btIt=dapResults.begin(); btIt != dapResults.end(); ++btIt) {
		    BaseType *bt = *btIt;
	        BESDEBUG("ugrid", "ugr4() - Adding variable "<< bt->name() << " to DAP structure " << dapResult->name() << endl);
			dapResult->add_var_nocopy(bt);
		}
	}

	*btpp = dapResult;

	BESDEBUG("ugrid", "ugr4() - Releasing maps and vectors..." << endl);
	for (mit = meshToRangeVarsMap->begin(); mit != meshToRangeVarsMap->end(); ++mit) {
	    vector<MeshDataVariable *> *requestedRangeVarsForMesh = mit->second;
        for(rvit=requestedRangeVarsForMesh->begin(); rvit!=requestedRangeVarsForMesh->end(); rvit++){
            MeshDataVariable *mdv = *rvit;
            delete mdv;
        }
		delete requestedRangeVarsForMesh;
	}
    delete meshToRangeVarsMap;





	BESDEBUG("ugrid", "ugr4() - END" << endl);

	return;
}


} // namespace ugrid_restrict
