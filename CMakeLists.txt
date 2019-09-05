cmake_minimum_required( VERSION 2.8.2 )

project( FLEXBAR )

set( SEQAN_APP_VERSION "3.5.0" )

include_directories( ${FLEXBAR_SOURCE_DIR}/include )
# link_directories( ${FLEXBAR_SOURCE_DIR}/lib )
#
include_directories(/home/linuxbrew/.linuxbrew/Cellar/tbb/2019_U8/include)
link_directories(/home/linuxbrew/.linuxbrew/Cellar/tbb/2019_U8/lib)

# file( MAKE_DIRECTORY build )

set( EXECUTABLE_OUTPUT_PATH ${FLEXBAR_BINARY_DIR} )

add_subdirectory( src )

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()
