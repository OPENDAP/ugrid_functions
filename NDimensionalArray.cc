// -*- mode: c++; c-basic-offset:4 -*-

// This file is part of libdap, A C++ implementation of the OPeNDAP Data
// Access Protocol.

// Copyright (c) 2002,2003,2011,2012 OPeNDAP, Inc.
// Authors: Nathan Potter <ndp@opendap.org>
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



#include "NDimensionalArray.h"
#include "util.h"

namespace libdap {



NDimensionalArray::NDimensionalArray()
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {

    string msg = "NDimArray::NDimArray() - INTERNAL_ERROR: This is the private constructor and should never be used";
    BESDEBUG(NDimensionalArray_debug_key, msg << endl);
    throw libdap::InternalErr(__FILE__, __LINE__, msg);
}


NDimensionalArray::NDimensionalArray(libdap::Array *a)
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {

    _shape = new vector<unsigned int>(a->dimensions(true), (unsigned int)1);
    _totalValueCount = computeConstrainedShape(a, _shape);
    _dapType = a->var()->type();

    allocateStorage(_totalValueCount, _dapType);
}


NDimensionalArray::NDimensionalArray(std::vector<unsigned int> *shape, libdap::Type dapType)
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {

    _shape = new vector<unsigned int>(*shape);
    _totalValueCount = computeArraySizeFromShapeVector(_shape);

    allocateStorage(_totalValueCount, _dapType);

}


NDimensionalArray::~NDimensionalArray() {
    delete (char *) _storage;
    delete _shape;
}

/**
 * Returns a pointer to the underlying storage for the NDimensionalArray. Calling this function is
 * effectively telling the instance of NDimensionalArray to 'release' it's reference to the storage. While
 * the memory will not be deleted by this call, the instance of NDimensionalArray will remove it's internal reference to
 * the storage and thus when the NDimensionalArray goes out of scope, or is otherwise deleted the storage WILL NOT BE DELETED.
 * CALLING THIS METHOD MEANS THAT YOU ARE NOW RESPONSIBLE FOR FREEING THE MEMORY REFERENCED BY THE RETURNED POINTER.
 */
void *NDimensionalArray::relinquishStorage(){
    void *s = _storage;
    _storage = 0;
    return s;
}



/**
 * Computes and returns (via the returned value parameter 'shape') the constrained shape of the libdap::Array 'a'.
 * Returns the total number of elements in constrained shape.
 */
long NDimensionalArray::computeConstrainedShape(libdap::Array *a, vector<unsigned int> *shape ){

    libdap::Array::Dim_iter dIt;
    unsigned int start;
    unsigned int stride;
    unsigned int stop;

    unsigned int dimSize = 1;
    int dimNum = 0;
    long totalSize = 1;

    for(dIt =a->dim_begin() ; dIt!=a->dim_end() ;dIt++){
        start  = a->dimension_start(dIt, true);
        stride = a->dimension_stride(dIt, true);
        stop   = a->dimension_stop(dIt, true);
        dimSize = 1 + ( (stop - start) / stride);
        (*shape)[dimNum++] = dimSize;
        totalSize *= dimSize;
    }

    return totalSize;
}

/**
 * Computes the total number of elements of the n-dimensional array described by the shape vector.
 */
long NDimensionalArray::computeArraySizeFromShapeVector(vector<unsigned int> *shape ){
    long totalSize = 1;

    for(int i; i<shape->size(); i++){
        totalSize *= (*shape)[i];
    }

    return totalSize;
}

/**
 * Computes the element index in the underlying one dimensional array for the passed location based on an
 * n-dimensional array described by the shape vector.
 */
long NDimensionalArray::getStorageIndex(vector<unsigned int> *shape, vector<unsigned int> *location){
    BESDEBUG(NDimensionalArray_debug_key, "NDimensionalArray::getStorageIndex() - BEGIN." << endl);
    long storageIndex = 0;


    if(location->size() != shape->size()){
        string msg = "getStorageIndex() - The supplied location vector does not match array shape.";
        BESDEBUG(NDimensionalArray_debug_key, msg << endl);
        throw Error(msg);
    }

    BESDEBUG(NDimensionalArray_debug_key, "NDimensionalArray::getStorageIndex() - Shape and location have the same number of elements." << endl);

    long dimIndex     = 0;
    long chunkSize = 1;

    for(dimIndex = shape->size()-1 ; dimIndex >=0 ; dimIndex--){
        BESDEBUG(NDimensionalArray_debug_key, "NDimensionalArray::getStorageIndex() - dimIndex=" << libdap::long_to_string(dimIndex)  << endl);

       if((*location)[dimIndex] >= (*shape)[dimIndex]){
            string msg = "NDimensionalArray::getStorageIndex() - The location vector references a value that does not match the array shape. ";
            msg += "location[" + libdap::long_to_string(dimIndex) + "]=";
            msg += libdap::long_to_string((*location)[dimIndex]) + " ";
            msg += "shape[" + libdap::long_to_string(dimIndex) + "]=";
            msg += libdap::long_to_string((*shape)[dimIndex]) + " ";
            BESDEBUG(NDimensionalArray_debug_key, msg << endl);
            throw Error(msg);
        }
        storageIndex += chunkSize * ((*location)[dimIndex]);
        chunkSize *= ((*shape)[dimIndex]);
    }

    BESDEBUG(NDimensionalArray_debug_key, "NDimensionalArray::getStorageIndex() - END." << endl);
    return storageIndex;
}

/**
 * Allocates internal storage for the NDimensionalArray
 */
void NDimensionalArray::allocateStorage(long numValues, Type dapType){

    switch (dapType) {
    case dods_byte_c:
        _sizeOfValue = sizeof(dods_byte);
        break;
    case dods_int16_c:
        _sizeOfValue = sizeof(dods_int16);
        break;
    case dods_uint16_c:
        _sizeOfValue = sizeof(dods_uint16);
        break;
    case dods_int32_c:
        _sizeOfValue = sizeof(dods_int32);
        break;
    case dods_uint32_c:
        _sizeOfValue = sizeof(dods_uint32);
        break;
    case dods_float32_c:
        _sizeOfValue = sizeof(dods_float32);
        break;
    case dods_float64_c:
        _sizeOfValue = sizeof(dods_float64);
        break;
    default:
        throw InternalErr(__FILE__, __LINE__,
                "Unknown DAP type encountered when constructing NDimensionalArray");
    }

    _storage = new char[numValues * _sizeOfValue];

}

/**
 * Verifies that the allocated storage for the NDimensioalArray has not been previously surrendered.
 */
void NDimensionalArray::confirmStorage(){
    if(_storage==0){
        string msg = "ERROR - NDimensionalArray storage has been relinquished. Instance is no longer viable for set/get operations.";
        BESDEBUG(NDimensionalArray_debug_key, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
}


/**
 * Verifies that the passed TypedapTypen is the same as the underlying type of the NDimensionalArray. If not, and Error is thrown.
 */
void NDimensionalArray::confirmType(Type dapType){
    if(_dapType != dapType){
        string msg = "NDimensionalArray::setValue() - Passed value does not match template array type. Expected "
                + libdap::type_name(_dapType) + " received "+ libdap::type_name(dapType);
        BESDEBUG(NDimensionalArray_debug_key, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
}

/**
 * Verifies that the passed value n is the same as the size of the last dimension. If not, and Error is thrown.
 */
void NDimensionalArray::confirmLastDimSize(unsigned int n){
    long elementCount = getLastDimensionElementCount();
    if(elementCount != n){
        string msg = "NDimensionalArray::setLastDimensionHyperSlab() - Passed valueCount does not match size of last dimension hyper-slab. ";
        msg += "Last dimension hyper-slab has " + libdap::long_to_string(elementCount) + " elements. ";
        msg += "Received a valueCount of  "+ libdap::long_to_string(n);
        BESDEBUG(NDimensionalArray_debug_key, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }

}




/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_byte NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_byte value){

    confirmStorage();
    confirmType(dods_byte_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_byte *_store = static_cast<dods_byte*>(_storage);
    dods_byte oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_int16 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_int16 value){

    confirmStorage();
    confirmType(dods_int16_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_int16 *_store = static_cast<dods_int16 *>(_storage);
    dods_int16 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_uint16 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_uint16 value){
    confirmStorage();
    confirmType(dods_uint16_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_uint16 *_store = static_cast<dods_uint16 *>(_storage);
    dods_uint16 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_int32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_int32 value){
    confirmStorage();
    confirmType(dods_int32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_int32 *_store = static_cast<dods_int32 *>(_storage);
    dods_int32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_uint32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_uint32 value){
    confirmStorage();
    confirmType(dods_uint32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_uint32 *_store = static_cast<dods_uint32 *>(_storage);
    dods_uint32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_float32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_float32 value){
    confirmStorage();
    confirmType(dods_float32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_float32 *_store = static_cast<dods_float32 *>(_storage);
    dods_float32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

/**
 * Sets the value of the array at the location specified by the passed location vector. The number
 * of elements in the location vector must be the the same as the number of dimensions in the NDimensionalArray
 * and the type of the passed value must match the underlying type of the NDimensionalArray
 */
dods_float64 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_float64 value){
    confirmStorage();
    confirmType(dods_float64_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_float64 *_store = static_cast<dods_float64 *>(_storage);
    dods_float64 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}


/**
 * The return value parameters slab and elementCount are used to return a pointer to the first element of the last dimension hyper-slab
 * and the number of elements in the hyper-slab. The passed the location vector, identifies the requested slab.The location vector must
 * have N-1 elements where N is the number of dimensions in the NDimensionalArray.
 */
void NDimensionalArray::getLastDimensionHyperSlab(std::vector<unsigned int> *location, void **slab, unsigned int *elementCount){
    confirmStorage();
    if(location->size()!=_shape->size()-1){
        string msg = "NDimensionalArray::getLastDimensionHyperSlab() - Passed location vector doesn't match array shape.";
        BESDEBUG(NDimensionalArray_debug_key, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
    vector<unsigned int> slabLocation(*location);

    slabLocation.push_back(0);
    unsigned int storageIndex = getStorageIndex(_shape, &slabLocation);

    *slab = &((char *)_storage)[storageIndex*_sizeOfValue];
    *elementCount = *(_shape->rbegin());

}



/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_byte *values, unsigned int valueCount){
    confirmType(dods_byte_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_byte));
}


/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int16 *values, unsigned int valueCount){
    confirmType(dods_int16_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_int16));
}


/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint16 *values, unsigned int valueCount){
    confirmType(dods_uint16_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_uint16));
}


/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int32 *values, unsigned int valueCount){
    confirmType(dods_int32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_int32));
}


/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint32 *values, unsigned int valueCount){
    confirmType(dods_uint32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_uint32));
}


/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float32 *values, unsigned int valueCount){
    confirmType(dods_float32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_float32));
}

/**
 * Sets all of the values in the last dimension hyper-hyper slab identified by the N-1 element location vector (where N is the
 * number of dimensions in the NDimensionalArray). The number of values passed in must match the size of the last dimension
 * hyper-slab, and the type of the values must match the underlying type of the NDimensionalArray.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float64 *values, unsigned int valueCount){
    confirmType(dods_float64_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_float64));
}



/**
 * This private method uses 'memcopy' to perform a byte by byte copy of the passed values array onto the last dimension
 * hyper-slab referenced by the N-1 element vector location.
 */
void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, void *values, unsigned int byteCount){
    confirmStorage();
    void *slab;
    unsigned int slabElementCount;

    getLastDimensionHyperSlab(location,&slab,&slabElementCount);
    memcpy ( slab, values, byteCount );

}

/**
 * Uses 'memset' to set ALL of the values in the NDimensionalArray to the passed char value.
 */
void NDimensionalArray::setAll(char val){
    confirmStorage();
    memset(_storage,val,_totalValueCount*_sizeOfValue);

}

/**
 * Returns the size, in elements, of the last dimension.
 */
long NDimensionalArray::getLastDimensionElementCount(){
    return *(_shape->rbegin());
}







} /* namespace libdap */
