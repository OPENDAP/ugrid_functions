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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// You can contact OPeNDAP, Inc. at PO Box 112, Saunderstown, RI. 02874-0112.

#include <iostream>

using std::endl;

#include "UgridFunctions.h"
#include "gf3.h"
#include "gf2.h"
#include "ce_functions.h"
#include "BESServerFunctionsList.h"
#include "BESDebug.h"

//namespace ugrid {

void UgridFunctions::initialize(const string &modname) {
	BESDEBUG( "UgridFunctions", "Initializing UgridFunctions:" << endl );

	BESServerFunctionsList::TheList()->add_function("ugr3", gf3::function_ugr3);
	BESServerFunctionsList::TheList()->add_function("ugr2", gf2::function_ugr2);

	BESDEBUG( "UgridFunctions", "Done initializing UgridFunctions" << endl );
}

void UgridFunctions::terminate(const string &modname) {
	BESDEBUG( "UgridFunctions", "Removing UgridFunctions Modules (this does nothing)." << endl );
}

/** @brief dumps information about this object
 *
 * Displays the pointer value of this instance
 *
 * @param strm C++ i/o stream to dump the information to
 */
void UgridFunctions::dump(ostream &strm) const {
	strm << BESIndent::LMarg << "UgridFunctions::dump - (" << (void *) this << ")" << endl;
}

extern "C" {
	BESAbstractModule *maker() {
		return new UgridFunctions;
	}
}




//} // namespace ugrid


