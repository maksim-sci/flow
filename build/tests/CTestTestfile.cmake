# CMake generated Testfile for 
# Source directory: D:/cpp/Sci/Flow/tests
# Build directory: D:/cpp/Sci/Flow/build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_electrode "D:/cpp/Sci/Flow/build/tests/test_electrode.exe")
set_tests_properties(test_electrode PROPERTIES  _BACKTRACE_TRIPLES "D:/cpp/Sci/Flow/tests/CMakeLists.txt;29;add_test;D:/cpp/Sci/Flow/tests/CMakeLists.txt;0;")
add_test(test_grid "D:/cpp/Sci/Flow/build/tests/test_grid.exe")
set_tests_properties(test_grid PROPERTIES  _BACKTRACE_TRIPLES "D:/cpp/Sci/Flow/tests/CMakeLists.txt;30;add_test;D:/cpp/Sci/Flow/tests/CMakeLists.txt;0;")
subdirs("../_deps/catch-build")
