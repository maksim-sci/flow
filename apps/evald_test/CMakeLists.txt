project(evald_test LANGUAGES CXX)

add_executable(evald_test main.cpp)

target_link_libraries(evald_test PRIVATE fmt::fmt flow_library)