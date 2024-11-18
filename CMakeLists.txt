#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(ice_spice)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------------
# Add option to debug
#
option(DEBUG_INTERPOLATING_FIELD "Debug setting" OFF)
if(DEBUG_INTERPOLATING_FIELD)
  add_definitions(-DDEBUG_INTERPOLATING_FIELD)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ICESPICE ICESPICE.cc ${sources} ${headers})
target_link_libraries(ICESPICE ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ICESPICE. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ICESPICE_SCRIPT
    vis.mac 
    gui.mac
    comsol_output_5N42_1x1x8in_x50_y50_z70_res1_2mm.mag
    comsol_output_3N42_1x1x16in_x50_y50_z70_res1_2mm.mag
    comsol_output_5N42_1x1x16in_x50_y50_z70_res1_2mm.mag
    comsol_output_6N42_1x1x16in_x50_y50_z70_res1_2mm.mag
    Plots.C
    Fold.C
  )

foreach(_script ${ICESPICE_SCRIPT})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# Copy the entire cad_files directory to the build directory
file(COPY ${PROJECT_SOURCE_DIR}/cad_files DESTINATION ${PROJECT_BINARY_DIR})

configure_file(${PROJECT_SOURCE_DIR}/207Bi/bi207.mac ${PROJECT_BINARY_DIR}/bi207.mac COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/207Pb/pb207.mac ${PROJECT_BINARY_DIR}/pb207.mac COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/TransmissionProbablilty/transmission_prob.mac ${PROJECT_BINARY_DIR}/transmission_prob.mac COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/BackScattering/backscattering.mac ${PROJECT_BINARY_DIR}/backscattering.mac COPYONLY)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ice_spice DEPENDS ICESPICE)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ICESPICE DESTINATION bin)