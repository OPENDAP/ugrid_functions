
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

#ifndef UGR4_H_
#define UGR4_H_

#include "BaseType.h"
#include "DDS.h"
#include "ServerFunction.h"

namespace ugrid {

void ugr4(int argc, libdap::BaseType * argv[], libdap::DDS &dds, libdap::BaseType **btpp) ;

/**
 * The UGR4 class encapsulates the function 'ugr4::ugr4'
 * along with additional meta-data regarding its use and applicability.
 */
class UGR4: public libdap::ServerFunction {

private:
#if 0
    string name;;
    string description;
    string usage;
    string role;
    string doumentationUrl;
    string version;
#endif

public:
    UGR4()
#if 0
:
        ServerFunction(),
        name("ugr4"),
        description("This function can subset the range variables of a two dimensional triangular mesh unstructured grid."),
        usage("ugr4(0, node_var [,node_var_2,...,node_var_n], 'relational query over range')"),
        role("http://services.opendap.org/dap4/server-side-function/unstructured_grids/ugrid_restrict"),
        doumentationUrl("http://docs.opendap.org/index.php/UGrid_Functions"),
        version("1.0")
#endif
{
#if 0
        setName(name);
        setDescriptionString(description);
        setUsageString(usage);
        setRole(role);
        setDocUrl(doumentationUrl);
        setFunction(ugrid::ugr4);
        setVersion(version);
#endif
#if 1
		setName("ugr4");
		setDescriptionString("This function can subset the range variables of a two dimensional triangular mesh unstructured grid.");
		setUsageString("ugr4(0, node_var [,node_var_2,...,node_var_n], 'relational query over range')");
		setRole("http://services.opendap.org/dap4/server-side-function/unstructured_grids/ugrid_restrict");
		setDocUrl("http://docs.opendap.org/index.php/UGrid_Functions");
		setFunction(ugrid::ugr4);
		setVersion("1.0");
#endif
    }
    virtual ~UGR4()
    {
    }

};


}// namespace ugrid_restrict


#endif /* UGR4_H_ */
