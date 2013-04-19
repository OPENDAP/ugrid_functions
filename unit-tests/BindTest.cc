// -*- mode: c++; c-basic-offset:4 -*-

// This file is part of libdap, A C++ implementation of the OPeNDAP Data
// Access Protocol.

// Copyright (c) 2005 OPeNDAP, Inc.
// Author: Nathan David Potter <ndp@opendap.org>
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
#include <cppunit/TextTestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cstdio>

#define DODS_DEBUG

#include "debug.h"


#include "gridfields/grid.h"
#include "gridfields/gridfield.h"
#include "gridfields/bind.h"
#include "gridfields/array.h"
#include "gridfields/restrict.h"
#include "gridfields/refrestrict.h"
#include "gridfields/arrayreader.h"
#include "gridfields/accumulate.h"

using namespace GF;

namespace ugrid {

class BindTest : public CppUnit::TestFixture {

private:

    Grid *makeGrid(int scale, string name) {
      CellArray *twocells;
      CellArray *onecells;
      CellArray *zerocells;
      Grid *grid;
      Node triangle[3];
      Node segment[2];
      Node node;

      bool wf;
      int i;
      twocells = new CellArray();
      for (i=0; i<scale/2; i++) {
        triangle[0] = i;
        triangle[1] = i+1;
        triangle[2] = i+2;
        twocells->addCellNodes(triangle, 3);
      }
      //twocells->print();
      //getchar();
      onecells = new CellArray();
      for (i=0; i<scale-1; i++) {
        segment[0] = i;
        segment[1] = i+1;
        onecells->addCellNodes(segment, 2);
      }
      //onecells->print();

      //getchar();
      grid = new Grid(name, 2);
      grid->setImplicit0Cells(scale);
      grid->setKCells(onecells, 1);
      grid->setKCells(twocells, 2);
      //grid->print(0);
      //getchar();
      return grid;
    }


    Array *makeFloatArray(int size,const char *name) {
      Array *arr;
      arr = new Array(name, FLOAT, size);
      float *data;
      arr->getData(data);
      int i;

      for (i=0; i<size; i++) {
          data[i] = 2*i-10;
      }
      return arr;
    }
    GridField *makeGridField(int size, string gridname,const char *datname, int k) {

      Grid *G;
      GridField *GF;
      Array *data;

      G = makeGrid(12, "A");
      k = 0;
      data = makeFloatArray(12, "x");

      GF = new GridField(G, k, data);
      //printf("Valid? %i\n", !notValid(GF));
      //GF->print();

      return GF;
    }


public:

    // Called once before everything gets tested
    BindTest() {
    //    DBG(cerr << " BindTest - Constructor" << endl);

    }

    // Called at the end of the test
    ~BindTest() {
    //    DBG(cerr << " BindTest - Destructor" << endl);
    }


    // Called before each test
    void setup() {
    //    DBG(cerr << " BindTest - setup()" << endl);
    }

    // Called after each test
    void tearDown() {
    //    DBG(cerr << " tearDown()" << endl);
    }

    CPPUNIT_TEST_SUITE( BindTest );

    CPPUNIT_TEST(bind_test);

    CPPUNIT_TEST_SUITE_END();


    void  bind_test() {
      bool verbose = true;

      try {
        GridField *GF;
        GridField *Result;

        GF = makeGridField(12, "A", "x", 0);
        Array *arr = new Array("io", FLOAT, 12);
        GF->Bind(0, arr);

        GridField *aGF = AccumulateOp::Accumulate(GF, 0, "result", "result+1", "0", 0);

        if (verbose) aGF->print(9);
        printf("restricting...\n");
        Result = RefRestrictOp::Restrict("x<4",0,GF);
        if (verbose) Result->print(0);
        Result = RefRestrictOp::Restrict("x>-4",0,Result);
        if (verbose) Result->print(10);

        FileArrayReader *ar = new FileArrayReader("bindtest.dat", 0);
        ar->setPatternAttribute("result");
        GridField *G = BindOp::Bind("io", FLOAT, ar, 0, Result);

        if (verbose) G->print();

        CPPUNIT_ASSERT(true);
      }
      catch (std::string &e) {
        cerr << "Error: " << e << endl;
        CPPUNIT_ASSERT(false);
      }
      catch (...) {
        cerr << "Unknown Error." << endl;
        CPPUNIT_ASSERT(false);
      }
    }

}; // BindTest

CPPUNIT_TEST_SUITE_REGISTRATION(BindTest);

} // namespace ugrid


int
main( int, char** )
{
    CppUnit::TextTestRunner runner;
    runner.addTest( CppUnit::TestFactoryRegistry::getRegistry().makeTest() );

    bool wasSuccessful = runner.run( "", false ) ;

    return wasSuccessful ? 0 : 1;
}


