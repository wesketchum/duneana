# ======================================================================
#  duneana main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(duneana VERSION 09.42.03.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_boost( )
find_ups_root()
find_ups_product( cetbuildtools )
find_ups_product( dunecore )
find_ups_product( dunereco )
find_ups_product( duneopdet )
find_ups_product( larreco )
find_ups_product( larsim )
find_ups_product( larevt )
find_ups_product( larpandora )
find_ups_product( systematicstools )
#if(DEFINED ENV{HEP_HPC_DIR} )
find_ups_product( clhep )
find_ups_product( hdf5 )
find_ups_product( hep_hpc )
find_ups_geant4( )


if(DEFINED ENV{HDF5_DIR} )
  string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
  if( ${BTYPE_UC} MATCHES "DEBUG" )
    cet_find_library( HDF5 NAMES hdf5_debug PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
  else()
    cet_find_library( HDF5 NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
  endif()
endif()


#endif()

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneana)
add_subdirectory(ups)

include(UseCPack)
