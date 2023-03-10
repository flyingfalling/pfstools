include_directories ("${PROJECT_BINARY_DIR}/"
"${PROJECT_SOURCE_DIR}/src/pfs" "${CMAKE_CURRENT_SOURCE_DIR}")

link_directories("${PROJECT_SOURCE_DIR}/src/pfs")

set(TRG pfshdrcalibrate)
add_executable(${TRG} ${TRG}.cpp  responses.cpp robertson02.cpp  mitsunaga99.cpp  mitsunaga99_numerical.cpp nrutil.cpp "${GETOPT_OBJECT}")			  
target_link_libraries(${TRG} pfs)
install (TARGETS ${TRG} DESTINATION bin)
install (FILES ${TRG}.1 DESTINATION ${MAN_DIR})

if( OpenCV_FOUND AND EXIF_FOUND )

	set(TRG pfsalign)
	add_executable(${TRG} ${TRG}.cpp  "${GETOPT_OBJECT}")			  
	target_link_libraries(${TRG} pfs ${OpenCV_LIBS}
	${EXIF_LIBRARIES})
	include_directories( ${EXIF_INCLUDE_DIRS} )
	install (TARGETS ${TRG} DESTINATION bin)
	install (FILES ${TRG}.1 DESTINATION ${MAN_DIR})

endif( OpenCV_FOUND AND EXIF_FOUND )


# Shell scripts
set(SCRIPTS dcraw2hdrgen jpeg2hdrgen pfsinme pfsplotresponse pfsinhdrgen) 

foreach(SCRIPT ${SCRIPTS}) 

	# Replace the tag with the path to bash
	file(READ ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT}.in file_content)
	string(REGEX REPLACE "@BASH_PATH@" "${BASH_EXECUTABLE}" file_content "${file_content}")
	string(REGEX REPLACE "@PERL_PATH@" "${PERL_EXECUTABLE}" file_content "${file_content}")
	file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}" "${file_content}")	

	install (FILES "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}" 
         PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE OWNER_WRITE WORLD_READ GROUP_READ OWNER_READ
         DESTINATION bin)	
		 
	install (FILES ${SCRIPT}.1 DESTINATION ${MAN_DIR})
		 
endforeach(SCRIPT) 

install (FILES pfs_split_exposures.py pfs_automerge  
	PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE OWNER_WRITE WORLD_READ GROUP_READ OWNER_READ
	DESTINATION bin)	

install (FILES pfs_automerge.1 DESTINATION ${MAN_DIR})


