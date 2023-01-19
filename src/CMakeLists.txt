if( NOT HAS_GETOPT )
  add_subdirectory (getopt)
endif( NOT HAS_GETOPT )

add_subdirectory (pfs) 
add_subdirectory (fileformat) 
add_subdirectory (filter)

if( ImageMagick_FOUND )
    add_subdirectory (hdrhtml)
endif( ImageMagick_FOUND )

if( Qt5Widgets_FOUND )
	add_subdirectory (pfsview)
else( Qt5Widgets_FOUND )
    # When pre-compiled pfsview is installed on Windows
    find_program( PFSVIEW pfsview )
    if( PFSVIEW )
    	message( "Found (precompiled) pfsview. pfsv script will be installed." )
        # Replace the tag with the path to bash
        file(READ ${CMAKE_CURRENT_SOURCE_DIR}/pfsview/pfsv.in file_content)
        string(REGEX REPLACE "@BASH_PATH@" "${BASH_EXECUTABLE}" file_content "${file_content}")
        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pfsview/pfsv" "${file_content}")	

        install (FILES "${CMAKE_CURRENT_BINARY_DIR}/pfsview/pfsv" 
            PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE OWNER_WRITE WORLD_READ GROUP_READ OWNER_READ
            DESTINATION bin)
    endif( PFSVIEW )
endif( Qt5Widgets_FOUND )

if( MATLAB_FOUND )
    add_subdirectory (matlab)
endif( MATLAB_FOUND )

if( OPENGL_FOUND AND GLUT_FOUND )
    add_subdirectory (pfsglview)
endif( OPENGL_FOUND AND GLUT_FOUND )

if( WITH_HDRVC )
	add_subdirectory (hdrvc)
endif( WITH_HDRVC )

if( MKOCTFILE )
	add_subdirectory (octave)
endif( MKOCTFILE )

add_subdirectory (tmo)
add_subdirectory (camera)
