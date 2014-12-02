
# 
# These macros are used for both the netcdf3 and netcdf4 tests.

AT_INIT([bes.conf besstandalone getdap])
# AT_COPYRIGHT([])

AT_TESTED([besstandalone])

AT_ARG_OPTION_ARG([generate g],
    [  -g arg, --generate=arg   Build the baseline file for test 'arg'],
    [if besstandalone -c $abs_builddir/bes.conf -i $abs_srcdir/$at_arg_generate \
     > $abs_srcdir/$at_arg_generate.baseline; then
         echo "Built baseline for $at_arg_generate"
     else
         echo "Could not generate baseline for $at_arg_generate"
     fi     
     exit],[])

AT_ARG_OPTION_ARG([data a],
    [  -a arg, --data=arg   Build the baseline file for test 'arg'],
    [if besstandalone -c $abs_builddir/bes.conf -i $abs_srcdir/$at_arg_data \
     | getdap -M - > $abs_srcdir/$at_arg_data.baseline; then
         echo "Built baseline for $at_arg_data"
     else
         echo "Could not generate baseline for $at_arg_data"
     fi     
     exit],[])

# Usage: _AT_TEST_*(<bescmd source>, <baseline file>)

m4_define([_AT_BESCMD_TEST],   
[AT_SETUP([BESCMD $1])
AT_KEYWORDS([bescmd])
AT_CHECK([besstandalone -c $abs_builddir/bes.conf -i $1 || true], [], [stdout], [stderr])
AT_CHECK([diff -b -B $2 stdout || diff -b -B $2 stderr], [], [ignore],[],[])
AT_CLEANUP])

m4_define([_AT_BESCMD_BINARYDATA_TEST],   
[AT_SETUP([BESCMD $1])
AT_KEYWORDS([bescmd])
AT_CHECK([besstandalone -c $abs_builddir/bes.conf -i $1 | getdap -M - || true], [], [stdout], [stderr])
AT_CHECK([diff -b -B $2 stdout || diff -b -B $2 stderr], [], [ignore],[],[])
AT_CLEANUP])

m4_define([AT_BESCMD_RESPONSE_TEST],
[_AT_BESCMD_TEST([$abs_srcdir/bescmd/$1], [$abs_srcdir/bescmd/$1.baseline])
])

m4_define([AT_BESCMD_BINARYDATA_RESPONSE_TEST],
[_AT_BESCMD_BINARYDATA_TEST([$abs_srcdir/bescmd/$1], [$abs_srcdir/bescmd/$1.baseline])
])
