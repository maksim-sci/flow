# Testing library
FetchContent_Declare(
        catch
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v2.13.6)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# Tests need to be added as executables first
add_executable(test_electrode testElectrode.cpp)
add_executable(test_grid testGridFromFile.cpp)

#fixing strange bug
if (MINGW)
  target_compile_options(test_electrode PRIVATE "-Wa,-mbig-obj")
  target_compile_options(test_grid PRIVATE "-Wa,-mbig-obj")
endif()

# I'm using C++17 in the test
target_compile_features(test_electrode PRIVATE cxx_std_17)
target_compile_features(test_grid PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(test_electrode PRIVATE vovka_library Catch2::Catch2)
target_link_libraries(test_grid PRIVATE vovka_library Catch2::Catch2)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME test_electrode COMMAND test_electrode) # Command can be a target
add_test(NAME test_grid COMMAND test_grid) # Command can be a target
