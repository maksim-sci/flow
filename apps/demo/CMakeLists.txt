project(demo LANGUAGES CXX)

add_executable(demo demo.cpp)

target_link_libraries(demo PRIVATE fmt::fmt flow_library)