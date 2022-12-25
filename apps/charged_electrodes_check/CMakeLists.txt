project(charged_electrodes_check LANGUAGES CXX)

add_executable(charged_electrodes_check main.cpp)

target_link_libraries(charged_electrodes_check PRIVATE fmt::fmt flow_library)