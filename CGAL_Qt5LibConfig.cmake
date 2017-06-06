set(WITH_CGAL_Qt5 "ON")

# The else condition of this code is never used in an installed
# version since it cannot happen there. Note also that for
# CMake<=2.8.11 (detected by the absence of CMP0024), the else()
# condition is never used.
if(NOT POLICY CMP0024 OR NOT CGAL_BUILDING_LIBS)
  if(NOT MSVC AND NOT CGAL_HEADER_ONLY)
    get_property(CGAL_Qt5_LIBRARY TARGET CGAL::CGAL_Qt5 PROPERTY LOCATION)
  else()
    set(CGAL_Qt5_LIBRARY "")
  endif()
else()
  # We are currently in a CGAL Build and CGALExports.cmake has not
  # necessarily been created yet. Just alias the targets. Also don't
  # access the LOCATION property here to set lib_LIBRARY, since those
  # targets are not imported and this is disallowed by CMP0026. Just
  # set it to the target name.
  if(TARGET CGAL_Qt5 AND NOT TARGET CGAL::CGAL_Qt5 AND NOT CGAL_HEADER_ONLY)
    add_library(CGAL::CGAL_Qt5 ALIAS CGAL_Qt5)
    set(CGAL_Qt5_LIBRARY CGAL::CGAL_Qt5)
  else()
    set(CGAL_Qt5_LIBRARY "")
  endif()
endif()


# 3RD_PARTY variables.
set(CGAL_Qt5_3RD_PARTY_INCLUDE_DIRS   "/usr/include")
set(CGAL_Qt5_3RD_PARTY_DEFINITIONS    "")
set(CGAL_Qt5_3RD_PARTY_LIBRARIES_DIRS "")
set(CGAL_Qt5_3RD_PARTY_LIBRARIES      "/usr/lib/x86_64-linux-gnu/libGLU.so;/usr/lib/x86_64-linux-gnu/libGL.so;/usr/lib/x86_64-linux-gnu/libSM.so;/usr/lib/x86_64-linux-gnu/libICE.so;/usr/lib/x86_64-linux-gnu/libX11.so;/usr/lib/x86_64-linux-gnu/libXext.so")
