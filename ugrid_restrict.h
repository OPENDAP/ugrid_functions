
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

#ifndef UGRID_RESTRICT_H_
#define UGRID_RESTRICT_H_

#include "BaseType.h"
#include "DDS.h"
#include "ServerFunction.h"

namespace ugrid_restrict {

void ugrid_restrict(int argc, libdap::BaseType * argv[], libdap::DDS &dds, libdap::BaseType **btpp) ;

/**
 * The UGrid3Function class encapsulates the echo arguments function 'gf3::function_ugr3'
 * along with additional meta-data regarding its use and applicability.
 */
class UGridRestrictFunction: public libdap::ServerFunction {
public:
	UGridRestrictFunction()
    {
		setName("ugr");
		setDescriptionString("This function can subset the range variables of a two dimensional triangular mesh unstructured grid.");
		setUsageString("ugr3(0, node_var [,node_var_2,...,node_var_n], 'relational query over range')");
		setRole("http://services.opendap.org/dap4/server-side-function/unstructured_grids/ugrid_restrict");
		setDocUrl("http://docs.opendap.org/index.php/UGrid_Functions");
		setFunction(ugrid_restrict::ugrid_restrict);
		setVersion("1.0");
    }
    virtual ~UGridRestrictFunction()
    {
    }

};


}// namespace ugrid_restrict


#endif /* UGRID_RESTRICT_H_ */
