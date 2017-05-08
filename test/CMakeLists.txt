cmake_minimum_required (VERSION 3.8)

set (HEADER_FILES common.hpp TestUtilities.hpp)

set (SOURCE_FILES colorsystem.cpp
                  spectrum.cpp
                  matrices.cpp
                  cielab.cpp
                  illuminants.cpp
                  aces2065.cpp
                  macbeth.cpp
                  units.cpp
    )

add_library (colortest_objs OBJECT ${SOURCE_FILES} ${HEADER_FILES})
target_compile_features (colortest_objs PUBLIC cxx_std_14)
target_include_directories (colortest_objs PUBLIC $<TARGET_PROPERTY:ColorSystem,INTERFACE_INCLUDE_DIRECTORIES>
                                                  ${COLORSYSTEM_SOURCE_DIR}/ext)
set_target_properties (colortest_objs PROPERTIES
                       COTIRE_PREFIX_HEADER_INIT common.hpp
                       COTIRE_PREFIX_HEADER_INCLUDE_PATH ${COLORSYSTEM_SOURCE_DIR}/ext)

if (COMMAND cotire)
    cotire (colortest_objs)
endif ()

add_executable (colortest main.cpp $<TARGET_OBJECTS:colortest_objs>)
target_link_libraries (colortest PRIVATE ColorSystem)
target_include_directories (colortest PRIVATE ${COLORSYSTEM_SOURCE_DIR}/ext)
target_compile_features (colortest PRIVATE cxx_std_14)

add_test (NAME colortest
          COMMAND colortest
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
