
# PROJECT_BINARY_DIR contains config.h and is thus needed
include_directories ("${PROJECT_SOURCE_DIR}/src/pfs" "${PROJECT_BINARY_DIR}")

if( CYGWIN )
	# Cygwin normally does not require getopt but matlab sourves will be compiled with an external matlab 
	# compiler, most probably Visual Studio
	include_directories ("${PROJECT_SOURCE_DIR}/src/getopt")
	set( SCRIPT_EXT ".bat" )
else( CYGWIN )
	set( SCRIPT_EXT "" )
endif( CYGWIN )

#include_directories ("${PROJECT_BINARY_DIR}/" "${PROJECT_SOURCE_DIR}/src/pfs" ${MATLAB_INCLUDE_DIR})
#link_directories("${PROJECT_SOURCE_DIR}/src/pfs")

execute_process(COMMAND ${MATLAB_ROOT}/bin/mexext${SCRIPT_EXT} OUTPUT_VARIABLE MEX_EXT OUTPUT_STRIP_TRAILING_WHITESPACE)

message( "Matlab: mexext: ${MEX_EXT}" )
set(SRCMEX pfsclose.cpp pfsopen.cpp pfspopen.cpp pfsput.cpp pfsget.cpp pfs_transform_colorspace.cpp pfspclose.cpp) 

set(SRCM Contents.m pfs_read_image.m pfs_read_luminance.m
pfs_read_rgb.m pfs_read_xyz.m pfs_shell.m pfs_test_shell.m
pfs_write_image.m pfs_write_luminance.m pfs_write_rgb.m
pfs_write_xyz.m pfsview.m pfs_fix_path.m pfs_exec.m) 

#cmake_policy(SET CMP0018 OLD)
#cmake_policy(SET CMP0042 OLD) This policy is not recognized in cmake version 2.8.12.2

SET( CMAKE_CXX_COMPILER ${MATLAB_MEX} )
SET( CMAKE_C_COMPILER   ${MATLAB_MEX} )
# Suffix and Prefix of the output target file
SET( CMAKE_SHARED_LIBRARY_SUFFIX ".${MEX_EXT}" )  # set suffix to .mexa64
SET( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix 
 
# Variables controlling the build-phrase
SET( CMAKE_CXX_FLAGS "-cxx -largeArrayDims" )
#SET( CMAKE_CXX_FLAGS "-cxx -largeArrayDims CXXFLAGS='$$CXXFLAGS -std=c++11'" )
#SET( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
 
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
set( INC_FLAGS "" )
foreach( IDIR ${inc_dirs} )
	set( INC_FLAGS "${INC_FLAGS} -I${IDIR}" )
endforeach( IDIR )
 
 
if( CYGWIN ) 
	# This is ugly but necessary to work around CMake limitations
	# The compiled objects need to be copied to ensure that the dependency check works
	# Linker expects the obejcts to be in the <OBJECT_DIR> so we need to have them there as well
	SET( CMAKE_CXX_COMPILE_OBJECT 
		"<CMAKE_CXX_COMPILER> <DEFINES> <FLAGS> ${INC_FLAGS} -outdir <OBJECT_DIR> -c <SOURCE>; cp <OBJECT_DIR>/$$(basename <SOURCE> .cpp ).obj <OBJECT>"
	)
else( CYGWIN )

  
  if( MEXGCC )
    set( GCC_COMMAND "GCC=${MEXGCC}" )
    MESSAGE( STATUS "Using gcc for MEX: ${MEXGCC}.\n" )
  else( MEXGCC )
    set( GCC_COMMAND "" )
  endif( MEXGCC )
  
  SET( CMAKE_CXX_COMPILE_OBJECT
		"<CMAKE_CXX_COMPILER> ${GCC_COMMAND} <DEFINES> <FLAGS> ${INC_FLAGS} -outdir <OBJECT_DIR> -c <SOURCE>; mv <OBJECT_DIR>/$$(basename <SOURCE> .cpp ).o <OBJECT>"    
#		"<CMAKE_CXX_COMPILER> <DEFINES> <FLAGS> ${INC_FLAGS} -outdir <OBJECT_DIR> -c <SOURCE>; mv <OBJECT_DIR>/$$(basename <SOURCE> .cpp ).o <OBJECT>"
#		"<CMAKE_CXX_COMPILER> <DEFINES> <FLAGS> ${INC_FLAGS} -outdir <OBJECT_DIR> -c <SOURCE>"
	)
endif( CYGWIN ) 

 
# Variables controlling the linking-phase
SET( CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS )    # remove -shared options. mex does not accept the "-shared" option
 
SET(
    CMAKE_CXX_CREATE_SHARED_LIBRARY
    "<CMAKE_CXX_COMPILER> ${GCC_COMMAND} -cxx <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> -output <TARGET> <OBJECTS> <LINK_LIBRARIES>"
#    "<CMAKE_CXX_COMPILER> GCC='/usr/bin/gcc-4.7' -cxx <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> -output <TARGET> <OBJECTS> <LINK_LIBRARIES> LDFLAGS='$LDFLAGS -static-libstdc++'"
)

#

#SET( CMAKE_SKIP_RPATH TRUE ) 
# Variables controlling the installation RPATH
SET( CMAKE_INSTALL_RPATH "\$ORIGIN" )

# CMake will reset RPATH at the installation phase, so we need to specify CMAKE_INSTALL_RPATH

SET( CMAKE_SKIP_RPATH TRUE ) 

MESSAGE( STATUS "mex.cmake is loaded.\n" )
 
# Since Matlab tends to use ancient version of stdc++, dynamic linkimg with pfs library 
# causes heaps of problems. It is better to statically add required files. 
SET( PFS_LIB_SRC ${PROJECT_SOURCE_DIR}/src/pfs/pfs.cpp ${PROJECT_SOURCE_DIR}/src/pfs/colorspace.cpp )
 
## To add an mex file, you can:
#  ADD_LIBRARY( myprogram SHARED myprogram.cxx )
## Make sure you specify it as SHARED, and make sure you use the cxx extension. 

add_definitions(-DMATLAB_MEX_FILE) 
foreach(SRC ${SRCMEX}) 
  get_filename_component(MEXNAME ${SRC} NAME_WE) 
  add_library(${MEXNAME} SHARED ${SRC} mex_utils.cpp ${PFS_LIB_SRC} ${MEXDEFFILE}) 
#  target_link_libraries(${MEXNAME} pfs ${MATLAB_LIBRARIES}) 
#  target_link_libraries(${MEXNAME} pfs ) 
#  set_target_properties(${MEXNAME} PROPERTIES PREFIX "" SUFFIX 
#".${MEX_EXT}") 

  set_property(TARGET ${MEXNAME} PROPERTY POSITION_INDEPENDENT_CODE OFF) # remove the -fPIC option. mex does not accept the "-fPIC" option

  install (TARGETS ${MEXNAME} DESTINATION ${MATLAB_DEST_DIR})
  install (FILES ${MEXNAME}.m DESTINATION ${MATLAB_DEST_DIR})

endforeach(SRC) 

foreach(SRC ${SRCM}) 
  install (FILES ${SRC} DESTINATION ${MATLAB_DEST_DIR})
endforeach(SRC) 
