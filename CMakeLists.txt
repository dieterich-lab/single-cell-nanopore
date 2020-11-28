cmake_minimum_required( VERSION 2.8.2 )

project( FLEXBAR )

set( SEQAN_APP_VERSION "3.5.0" )

include_directories( ${FLEXBAR_SOURCE_DIR}/seqan/include )
# link_directories( ${FLEXBAR_SOURCE_DIR}/lib )

# file( MAKE_DIRECTORY build )

set( EXECUTABLE_OUTPUT_PATH "../bin" )

add_subdirectory( src )

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()
