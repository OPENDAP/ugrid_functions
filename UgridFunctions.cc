// UgridFunctions.cc

// This file is part of bes, A C++ back-end server implementation framework
// for the OPeNDAP Data Access Protocol.

// Copyright (c) 2013 OPeNDAP, Inc.
// Author: James Gallagher <jgallagher@opendap.org>
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

#include <iostream>

using std::endl;

#include "UgridFunctions.h"
#include "ServerFunctionsList.h"
#include "BESDebug.h"
#include "ugr5.h"

static string getFunctionNames()
{
    vector<string> names;
    libdap::ServerFunctionsList::TheList()->getFunctionNames(&names);

    string msg;
    for (std::vector<string>::iterator it = names.begin(); it != names.end(); ++it) {
        if (!msg.empty()) msg += ", ";

        msg += *it;
    }
    return msg;
}
void UgridFunctions::initialize(const string &/*modname*/)
{
    BESDEBUG("UgridFunctions", "initialize() - BEGIN" << endl);
    BESDEBUG("UgridFunctions", "initialize() - function names: " << getFunctionNames()<< endl);

#if 0
    BESDEBUG("UgridFunctions", "initialize() - Adding gf3::UGridRestrictFunction_03()" << endl);
    libdap::ServerFunctionsList::TheList()->add_function(new gf3::UGridRestrictFunction_03());
    BESDEBUG("UgridFunctions", "initialize() - function names: " << getFunctionNames()<< endl);

    BESDEBUG("UgridFunctions", "initialize() - Adding ugrid_restrict::UGridRestrictFunction()" << endl);
    libdap::ServerFunctionsList::TheList()->add_function(new ugrid_restrict::UGridRestrictFunction());
    BESDEBUG("UgridFunctions", "initialize() - function names: " << getFunctionNames()<< endl);
    BESDEBUG("UgridFunctions", "initialize() - Adding UGR4 function..." << endl);
    ugrid::UGR4 *ugr4 = new ugrid::UGR4();
    libdap::ServerFunctionsList::TheList()->add_function(ugr4);
    BESDEBUG("UgridFunctions", "initialize() - function names: " << getFunctionNames()<< endl);
#endif

    BESDEBUG("UgridFunctions", "initialize() - Adding UGR5 function..." << endl);

    ugrid::UGR5 *ugr5 = new ugrid::UGR5();
    libdap::ServerFunctionsList::TheList()->add_function(ugr5);

    BESDEBUG("UgridFunctions", "initialize() - function names: " << getFunctionNames()<< endl);

    BESDEBUG("UgridFunctions", "initialize() - END" << endl);
}

void UgridFunctions::terminate(const string &/*modname*/)
{
    BESDEBUG("UgridFunctions", "Removing UgridFunctions Modules (this does nothing)." << endl);
}

/** @brief dumps information about this object
 *
 * Displays the pointer value of this instance
 *
 * @param strm C++ i/o stream to dump the information to
 */
void UgridFunctions::dump(ostream &strm) const
{
    strm << BESIndent::LMarg << "UgridFunctions::dump - (" << (void *) this << ")" << endl;
}

extern "C" {
BESAbstractModule *maker()
{
    return new UgridFunctions;
}
}
