cmake_minimum_required(VERSION 3.11)

project(
        flow
        VERSION 0.11
        DESCRIPTION "ex kmmk"
        LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    include(CTest)
endif()

  set(CMAKE_CXX_STANDARD 20)

include(FetchContent)


link_libraries(stdc++fs)


# Formatting library
FetchContent_Declare(
        fmtlib
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 5.3.0)
FetchContent_MakeAvailable(fmtlib)
# Adds fmt::fmt

# The compiled library code is here
add_subdirectory(src)


# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
# if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
#         AND BUILD_TESTING)
#     add_subdirectory(tests)
# endif()



#================================
#================================
#================================
#================================
#================================
#================================
#TESTING

include(CTest) 


enable_testing()

