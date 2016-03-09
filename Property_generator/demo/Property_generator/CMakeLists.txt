###############################################################################

# Add a new module to cmake to find Qwt. 
# If Qwt is not being found, try modifying the search directories
# defined in this module.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

cmake_minimum_required(VERSION 2.8)

###############################################################################

# The project. 
project( delaunay_visualiser )


#include_directories(${QWT_INCLUDE_DIRS})
#set(LIBS ${LIBS} ${QWT_LIBRARIES})

###############################################################################

FIND_PACKAGE(Qt4                 REQUIRED)
FIND_PACKAGE(QWT                 REQUIRED)
FIND_PACKAGE(CGAL COMPONENTS Qt4 REQUIRED)
FIND_PACKAGE(OPENGL              REQUIRED)
	
	
	
	
if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND QWT_FOUND)

	INCLUDE(${CGAL_USE_FILE})	
    INCLUDE(${QT_USE_FILE})	

    INCLUDE_DIRECTORIES(${QWT_INCLUDE_DIR} ${QT_QTOPENGL_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
	INCLUDE_DIRECTORIES( ../../include )	
	INCLUDE_DIRECTORIES( ./include )
	
	SET( QT_USE_QTXML     TRUE )
	SET( QT_USE_QTMAIN    TRUE )
	SET( QT_USE_QTSCRIPT  TRUE )
	SET( QT_USE_QTOPENGL  TRUE )
	SET( CMAKE_CXX_FLAGS "-Wall -O3 -g")

	SET( SOURCE_DIR  ./src/     )
	SET( INCLUDE_DIR ./include/ )
	SET( UI_DIR      ./ui/      )
	

	# WALK SIMULATE #
	
    SET( walk_simulate_SOURCES ${SOURCE_DIR}/simulate.cpp )
	
    SET( walk_simulate_HEADERS 
		 ${INCLUDE_DIR}/simulate.h
	     ${INCLUDE_DIR}/walk.h )

	# WALK VISUALISATION #

    SET( delaunay_visualiser_SOURCES 
		 ${SOURCE_DIR}/main.cpp  
		 ${SOURCE_DIR}/visualiser.cpp 
		 ${SOURCE_DIR}/statistics_plotter.cpp 
		 ${SOURCE_DIR}/point_generator.cpp )

    SET( delaunay_visualiser_HEADERS 
 		 ${INCLUDE_DIR}/walk.h 
		 ${INCLUDE_DIR}/vertex_walk.h		
		 ${INCLUDE_DIR}/triangulation_viewer.h	
		 ${INCLUDE_DIR}/walk_traits.h 
		 ${INCLUDE_DIR}/quad_tree.h )
	
	SET( delaunay_visualiser_MOC_HEADERS 
		 ${INCLUDE_DIR}/visualiser.h 
		 ${INCLUDE_DIR}/statistics_plotter.h 
		 ${INCLUDE_DIR}/point_generator.h )
	
	SET( delaunay_visualiser_UIS 
		 ${UI_DIR}/visualiser.ui 
		 ${UI_DIR}/point_generator.ui )
		
		
    QT4_WRAP_UI( delaunay_visualiser_GEN_UIS  ${delaunay_visualiser_UIS} )
		
	# Note Patch to stop MOC tripping up. 
	QT4_WRAP_CPP(delaunay_visualiser_GEN_MOC_HEADERS 
		${delaunay_visualiser_MOC_HEADERS} 
		OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED )	    

	ADD_DEFINITIONS(${QT_DEFINITIONS})	

    set( MACOSX_RESOURCE_FILES ${CMAKE_SOURCE_DIR}/resources/walk_visualiser.icns ) 


	IF(APPLE)
	  # set how it shows up in the Info.plist file
	  SET(MACOSX_BUNDLE_ICON_FILE walk_visualiser.icns) 
	  # set where in the bundle to put the icns file
	  SET_SOURCE_FILES_PROPERTIES(./resources/walk_visualiser.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
	  # include the icns file in the target
	  SET(delaunay_visualiser_SOURCES ${delaunay_visualiser_SOURCES} ./resources/walk_visualiser.icns)
	ENDIF(APPLE)


	# Walk simulation terminal application
    ADD_EXECUTABLE( walk_simulate 
		${walk_simulate_SOURCES} 
		${walk_simulate_HEADERS} )

	# Walk visualisation application
    ADD_EXECUTABLE( "Delaunay Visualiser" MACOSX_BUNDLE
		${delaunay_visualiser_SOURCES}  
		${delaunay_visualiser_HEADERS}  
		${delaunay_visualiser_GEN_MOC_HEADERS} 
		${delaunay_visualiser_GEN_UIS} )
        
    # Link with libraries
    TARGET_LINK_LIBRARIES( "Delaunay Visualiser" ${QWT_LIBRARIES} )
    #target_link_libraries( "Delaunay Visualiser" ${QT_LIBRARIES} )
    TARGET_LINK_LIBRARIES( "Delaunay Visualiser" ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})
	  
	SET(plugin_dest_dir bin)
	SET(qtconf_dest_dir bin)
	SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/Delaunay\ Visualiser")
	IF(APPLE)
	  SET(plugin_dest_dir Delaunay\ Visualiser.app/Contents/MacOS)
	  SET(qtconf_dest_dir Delaunay\ Visualiser.app/Contents/Resources)
	  SET(APPS "\${CMAKE_INSTALL_PREFIX}/Delaunay\ Visualiser.app")
	ENDIF(APPLE)	  
	  
	INSTALL(TARGETS "Delaunay Visualiser" 
	    BUNDLE DESTINATION . COMPONENT Runtime
	    RUNTIME DESTINATION bin COMPONENT Runtime
	    )

	INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/imageformats" DESTINATION ${plugin_dest_dir}/plugins COMPONENT Runtime)

	INSTALL(CODE "
	    file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
	    " COMPONENT Runtime)

	SET(DIRS 
		${QT_LIBRARY_DIRS}
		${CGAL_USE_FILE}		 
		${CGAL_LIBRARIES}
		${CGAL_3RD_PARTY_LIBRARIES}
		${QWT_INCLUDE_DIR}
		${QWT_LIBS}
		# This seems to be needed on my machine... TODO: fix it.
		/Users/rhemsley/qwt-6.0.1/lib/
		/usr/local/qwt-6.0.1/lib/
	 	${QT_QTOPENGL_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
	
	INSTALL(CODE "
	    file(GLOB_RECURSE QTPLUGINS
	      \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
	    include(BundleUtilities)
	    fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${DIRS}\")
	    " COMPONENT Runtime)

	# To Create a package, one can run "cpack -G DragNDrop CPackConfig.cmake" on Mac OS X
	# where CPackConfig.cmake is created by including CPack
	# And then there's ways to customize this as well
	SET(CPACK_BINARY_DRAGNDROP ON)
	INCLUDE(CPack)

	IF(APPLE)
		MESSAGE(STATUS "NOTE: run 'cpack -G DragNDrop CPackConfig.cmake' to build an installer for mac.")
	ENDIF(APPLE)

else()

    MESSAGE(STATUS "This program requires the CGAL, QT and Qwt")
  
endif()
###############################################################################