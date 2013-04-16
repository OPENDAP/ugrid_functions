/*
 * NDimensionalArray.cc
 *
 *  Created on: Apr 10, 2013
 *      Author: ndp
 */

#include "NDimensionalArray.h"
#include "util.h"

namespace libdap {



NDimensionalArray::NDimensionalArray()
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {

    debugKey = NDimensionalArray_debug_key;
    string msg = "NDimArray::NDimArray() - INTERNAL_ERROR: This is the private constructor and should never be used";
    BESDEBUG(debugKey, msg << endl);
    throw libdap::InternalErr(__FILE__, __LINE__, msg);
}


NDimensionalArray::NDimensionalArray(libdap::Array *a)
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {
    debugKey = NDimensionalArray_debug_key;

    _shape = new vector<unsigned int>(a->dimensions(true), (unsigned int)1);
    _totalValueCount = computeConstrainedShape(a, _shape);
    _dapType = a->var()->type();

    allocateStorage(_totalValueCount, _dapType);
}


NDimensionalArray::NDimensionalArray(std::vector<unsigned int> *shape, libdap::Type dapType)
    :_storage(0),_totalValueCount(0),_shape(0),_sizeOfValue(0),_dapType(dods_null_c) {
    debugKey = NDimensionalArray_debug_key;

    _shape = new vector<unsigned int>(*shape);
    _totalValueCount = computeArraySizeFromShapeVector(_shape);

    allocateStorage(_totalValueCount, _dapType);

}


NDimensionalArray::~NDimensionalArray() {
    delete _storage;
}

void *NDimensionalArray::relinquishStorage(){
    void *s = _storage;
    _storage = 0;
    return s;
}




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

long NDimensionalArray::computeArraySizeFromShapeVector(vector<unsigned int> *shape ){
    long totalSize = 1;

    for(int i; i<shape->size(); i++){
        totalSize *= (*shape)[i];
    }

    return totalSize;
}

long NDimensionalArray::getStorageIndex(vector<unsigned int> *shape, vector<unsigned int> *location){
    //cout << "NDimensionalArray::getStorageIndex() - BEGIN." << endl;
    long storageIndex = 0;

    string debugKey = "ugrid";

    if(location->size() != shape->size()){
        string msg = "getStorageIndex() - The supplied location vector does not match array shape.";
        BESDEBUG(debugKey, msg << endl);
        throw Error(msg);
    }

    //cout << "NDimensionalArray::getStorageIndex() - Shape and location have the same number of elements." << endl;

    long dimIndex     = 0;
    long chunkSize = 1;

    for(dimIndex = shape->size()-1 ; dimIndex >=0 ; dimIndex--){
       //cout << "NDimensionalArray::getStorageIndex() - dimIndex="<< libdap::long_to_string(dimIndex) << endl;
       if((*location)[dimIndex] >= (*shape)[dimIndex]){
            string msg = "NDimensionalArray::getStorageIndex() - The location vector references a value that does not match the array shape. ";
            msg += "location[" + libdap::long_to_string(dimIndex) + "]=";
            msg += libdap::long_to_string((*location)[dimIndex]) + " ";
            msg += "shape[" + libdap::long_to_string(dimIndex) + "]=";
            msg += libdap::long_to_string((*shape)[dimIndex]) + " ";
            BESDEBUG(debugKey, msg << endl);
            throw Error(msg);
        }
        storageIndex += chunkSize * ((*location)[dimIndex]);
        chunkSize *= ((*shape)[dimIndex]);
    }

    //cout << "NDimensionalArray::getStorageIndex() - END." << endl;
    return storageIndex;
}

void NDimensionalArray::allocateStorage(long numValues, Type dapType){

    switch (dapType) {
    case dods_byte_c:
    {
        _sizeOfValue = sizeof(dods_byte);
        _storage = new dods_byte[_totalValueCount];
        break;
    }
    case dods_int16_c:
    {
        _sizeOfValue = sizeof(dods_int16);
        _storage = new dods_int16[_totalValueCount];
        break;
    }
    case dods_uint16_c:
    {
        _sizeOfValue = sizeof(dods_uint16);
        _storage = new dods_uint16[_totalValueCount];
        break;
    }
    case dods_int32_c:
    {
        _sizeOfValue = sizeof(dods_int32);
        _storage = new dods_int32[_totalValueCount];
        break;
    }
    case dods_uint32_c:
    {
        _sizeOfValue = sizeof(dods_uint32);
        _storage = new dods_uint32[_totalValueCount];
        break;
    }
    case dods_float32_c:
    {
        _sizeOfValue = sizeof(dods_float32);
        _storage = new dods_float32[_totalValueCount];
        break;
    }
    case dods_float64_c:
    {
        _sizeOfValue = sizeof(dods_float64);
        _storage = new dods_float64[_totalValueCount];
        break;
    }
    default:
        throw InternalErr(__FILE__, __LINE__,
                "Unknown DAP type encountered when constructing NDimensionalArray");
    }

}

void NDimensionalArray::confirmStorage(){
    if(_storage==0){
        string msg = "ERROR - NDimensionalArray storage has been relinquished. Instance is no longer viable for set/get operations.";
        BESDEBUG(debugKey, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
}


void NDimensionalArray::confirmType(Type dapType){
    if(_dapType != dapType){
        string msg = "NDimensionalArray::setValue() - Passed value does not match template array type. Expected "
                + libdap::type_name(_dapType) + " received "+ libdap::type_name(dapType);
        BESDEBUG(debugKey, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
}
void NDimensionalArray::confirmLastDimSize(unsigned int n){
    unsigned int elementCount = *(_shape->rbegin());
    if(elementCount != n){
        string msg = "NDimensionalArray::setLastDimensionHyperSlab() - Passed valueCount does not match size of last dimension hyper-slab. ";
        msg += "Last dimension hyper-slab has " + libdap::long_to_string(elementCount) + " elements. ";
        msg += "Received a valueCount of  "+ libdap::long_to_string(n);
        BESDEBUG(debugKey, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }

}





dods_byte NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_byte value){

    confirmStorage();
    confirmType(dods_byte_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_byte *_store = static_cast<dods_byte*>(_storage);
    dods_byte oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_int16 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_int16 value){

    confirmStorage();
    confirmType(dods_int16_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_int16 *_store = static_cast<dods_int16 *>(_storage);
    dods_int16 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_uint16 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_uint16 value){
    confirmStorage();
    confirmType(dods_uint16_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_uint16 *_store = static_cast<dods_uint16 *>(_storage);
    dods_uint16 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_int32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_int32 value){
    confirmStorage();
    confirmType(dods_int32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_int32 *_store = static_cast<dods_int32 *>(_storage);
    dods_int32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_uint32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_uint32 value){
    confirmStorage();
    confirmType(dods_uint32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_uint32 *_store = static_cast<dods_uint32 *>(_storage);
    dods_uint32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_float32 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_float32 value){
    confirmStorage();
    confirmType(dods_float32_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_float32 *_store = static_cast<dods_float32 *>(_storage);
    dods_float32 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}

dods_float64 NDimensionalArray::setValue(std::vector<unsigned int> *location, dods_float64 value){
    confirmStorage();
    confirmType(dods_float64_c);

    unsigned int storageIndex = getStorageIndex(_shape,location);
    dods_float64 *_store = static_cast<dods_float64 *>(_storage);
    dods_float64 oldValue = _store[storageIndex];
    _store[storageIndex] = value;
    return oldValue;
}


void NDimensionalArray::getLastDimensionHyperSlab(std::vector<unsigned int> *location, void **slab, unsigned int *elementCount){
    confirmStorage();
    if(location->size()!=_shape->size()-1){
        string msg = "NDimensionalArray::getLastDimensionHyperSlab() - Passed location vector doesn't match array shape.";
        BESDEBUG(debugKey, msg << endl);
        throw InternalErr(__FILE__, __LINE__, msg);
    }
    vector<unsigned int> slabLocation(*location);

    slabLocation.push_back(0);
    unsigned int storageIndex = getStorageIndex(_shape, &slabLocation);

    *slab = &((char *)_storage)[storageIndex*_sizeOfValue];
    *elementCount = *(_shape->rbegin());

}



void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_byte *values, unsigned int valueCount){
    confirmType(dods_byte_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_byte));
}


void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int16 *values, unsigned int valueCount){
    confirmType(dods_int16_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_int16));
}


void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint16 *values, unsigned int valueCount){
    confirmType(dods_uint16_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_uint16));
}


void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int32 *values, unsigned int valueCount){
    confirmType(dods_int32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_int32));
}


void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint32 *values, unsigned int valueCount){
    confirmType(dods_uint32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_uint32));
}


void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float32 *values, unsigned int valueCount){
    confirmType(dods_float32_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_float32));
}

void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float64 *values, unsigned int valueCount){
    confirmType(dods_float64_c);
    confirmLastDimSize(valueCount);
    setLastDimensionHyperSlab(location, (void *) values, valueCount*sizeof(dods_float64));
}



void NDimensionalArray::setLastDimensionHyperSlab(std::vector<unsigned int> *location, void *values, unsigned int byteCount){
    confirmStorage();
    void *slab;
    unsigned int slabElementCount;

    getLastDimensionHyperSlab(location,&slab,&slabElementCount);
    memcpy ( slab, values, byteCount );

}






} /* namespace libdap */
