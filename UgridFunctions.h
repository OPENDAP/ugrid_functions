// UgridFunctions.h

// This file is part of bes, A C++ back-end server implementation framework
// for the OPeNDAP Data Access Protocol.

// Copyright (c) 2013 OPeNDAP, Inc.
// Author: Nathan Potter <ndp@opendap.org>
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


#ifndef UGRIDFUNCTIONS_H_
#define UGRIDFUNCTIONS_H_

//namespace ugrid {

#include "BESAbstractModule.h"
#include "AbstractFunction.h"
#include "gf3.h"

class UgridFunctions: public BESAbstractModule {
public:
	UgridFunctions()
    {
    }
    virtual ~UgridFunctions()
    {
    }
    virtual void initialize(const string &modname);
    virtual void terminate(const string &modname);

    virtual void dump(ostream &strm) const;
};


/**
 * The UGrid3Function class encapsulates the echo arguments function 'gf3::function_ugr3'
 * along with additional meta-data regarding its use and applicability.
 */
class UGridRestrict3Function: public libdap::AbstractFunction {
public:
	UGridRestrict3Function()
    {
		setName("ugr3");
		setDescriptionString("This function can subset the node data of a two dimensional triangular mesh unstructured grid.");
		setUsageString("ugr3(0, node_var [,node_var_2,...,node_var_n], 'relational query over range')");
		setRole("http://services.opendap.org/dap4/server-side-function/unstructured_grids/ugrid_restrict");
		setDocUrl("http://docs.opendap.org/index.php/UGrid_Functions");
		setFunction(gf3::function_ugr3);
		setVersion("1.0");
    }
    virtual ~UGridRestrict3Function()
    {
    }

};




//} /* namespace ugrid */
#endif /* UGRIDFUNCTIONS_H_ */
