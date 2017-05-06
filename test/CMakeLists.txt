cmake_minimum_required (VERSION 3.8)

add_executable (colortest colorsystem.cpp)
target_compile_features (colortest PRIVATE cxx_std_14)
target_link_libraries (colortest PRIVATE ColorSystem)

add_test (NAME colortest
          COMMAND colortest
          WORKING DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
