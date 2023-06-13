project(perftests LANGUAGES CXX)

add_executable(perftests runtests.cpp)

target_link_libraries(perftests PRIVATE fmt::fmt flow_library)

add_compile_options(-Wmultichar)

add_test(NAME insert_check COMMAND perftests insert_check)
add_test(NAME check_iteration COMMAND perftests check_iteration)
add_test(NAME check_lookup COMMAND perftests check_lookup)


#add_test(NAME test_must_crush COMMAND perftests mc)