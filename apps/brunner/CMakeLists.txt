project(basic_runner LANGUAGES CXX)

add_executable(basic_runner main.cpp)

target_link_libraries(basic_runner PRIVATE fmt::fmt flow_library)