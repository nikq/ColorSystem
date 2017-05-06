cmake_minimum_required (VERSION 3.8)

add_executable (colortest main.cpp colorsystem.cpp)
target_compile_features (colortest PRIVATE cxx_std_14)
target_link_libraries (colortest PRIVATE ColorSystem)
target_include_directories (colortest PRIVATE ${COLORSYSTEM_SOURCE_DIR}/ext)
set_target_properties (colortest PROPERTIES
                       COTIRE_PREFIX_HEADER_INIT common.hpp
                       COTIRE_PREFIX_HEADER_INCLUDE_PATH ${COLORSYSTEM_SOURCE_DIR}/ext)

if (COMMAND cotire)
    cotire (colortest)
endif ()

add_test (NAME colortest
          COMMAND colortest
          WORKING DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
