/*
 * NDimensionalArray.h
 *
 *  Created on: Apr 10, 2013
 *      Author: ndp
 */

#ifndef NDIMENSIONALARRAY_H_
#define NDIMENSIONALARRAY_H_

#include "BESDebug.h"

#include "Array.h"



namespace libdap {

static string NDimensionalArray_debug_key = "ugrid";

class NDimensionalArray {
private:

    string debugKey;

    NDimensionalArray();

    libdap::Type _dapType;

    std::vector<unsigned int> *_shape;
    long _totalValueCount; // Number of elements
    unsigned int _sizeOfValue;
    void *_storage;

    void allocateStorage(long numValues, libdap::Type dapType);
    void confirmStorage();
    void confirmType(Type dapType);
    void confirmLastDimSize(unsigned int n);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, void *values, unsigned int byteCount);

public:

    NDimensionalArray(libdap::Array *arrayTemplate);
    NDimensionalArray(std::vector<unsigned int> *shape, libdap::Type dapType);

    virtual ~NDimensionalArray();

    dods_byte    setValue(std::vector<unsigned int> *location, dods_byte    value);
    dods_int16   setValue(std::vector<unsigned int> *location, dods_int16   value);
    dods_uint16  setValue(std::vector<unsigned int> *location, dods_uint16  value);
    dods_int32   setValue(std::vector<unsigned int> *location, dods_int32   value);
    dods_uint32  setValue(std::vector<unsigned int> *location, dods_uint32  value);
    dods_float32 setValue(std::vector<unsigned int> *location, dods_float32 value);
    dods_float64 setValue(std::vector<unsigned int> *location, dods_float64 value);

    static long computeConstrainedShape(libdap::Array *a, vector<unsigned int> *shape );
    static long computeArraySizeFromShapeVector(vector<unsigned int> *shape );
    static long getStorageIndex(vector<unsigned int> *shape, vector<unsigned int> *location);

    long elementCount() { return _totalValueCount; }
    unsigned int sizeOfElement() { return _sizeOfValue; }

    void *relinquishStorage();

    void *getStorage() { return _storage; }

    Type getTypeTemplate(){ return _dapType; }

    void getLastDimensionHyperSlab(std::vector<unsigned int> *location,void **slab, unsigned int *elementCount);

    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_byte    *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int16   *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint16  *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_int32   *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_uint32  *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float32 *values, unsigned int numVal);
    void setLastDimensionHyperSlab(std::vector<unsigned int> *location, dods_float64 *values, unsigned int numVal);


}; //NdimensionalArray

} /* namespace libdap */
#endif /* NDIMENSIONALARRAY_H_ */
