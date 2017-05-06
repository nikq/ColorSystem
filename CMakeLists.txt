cmake_minimum_required (VERSION 3.8)

project (COLORSYSTEM
         LANGUAGES CXX C
         VERSION 0.0.1)

include (GNUInstallDirs)
add_library (ColorSystem INTERFACE)
target_include_directories (ColorSystem INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                            $<INSTALL_INTERFACE:include>)

enable_testing ()

add_subdirectory (test)

install (DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
