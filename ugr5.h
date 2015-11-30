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

#ifndef UGR5_H_
#define UGR5_H_

#include "BaseType.h"
#include "DDS.h"
#include "ServerFunction.h"

namespace ugrid {

void ugr5(int argc, libdap::BaseType * argv[], libdap::DDS &dds, libdap::BaseType **btpp);

/**
 * The UGR4 class encapsulates the function 'ugr4::ugr4'
 * along with additional meta-data regarding its use and applicability.
 */
class UGR5: public libdap::ServerFunction {

private:

public:
    UGR5()
    {
        setName("ugr5");
        setDescriptionString(
            "This function can subset the range variables of a two dimensional triangular mesh unstructured grid.");
        setUsageString("ugr5(0, node_var [,node_var_2,...,node_var_n], 'relational query over domain')");
        setRole("http://services.opendap.org/dap4/server-side-function/unstructured_grids/ugrid_restrict");
        setDocUrl("http://docs.opendap.org/index.php/UGrid_Functions");
        setFunction(ugrid::ugr5);
        setVersion("1.0");
    }
    virtual ~UGR5()
    {
    }

};

} // namespace ugrid_restrict

#endif /* UGR5_H_ */
