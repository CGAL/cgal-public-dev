set(WITH_CGAL "ON")

# The else condition of this code is never used in an installed
# version since it cannot happen there. Note also that for
# CMake<=2.8.11 (detected by the absence of CMP0024), the else()
# condition is never used.
if(NOT POLICY CMP0024 OR NOT CGAL_BUILDING_LIBS)
  if(NOT MSVC AND NOT CGAL_HEADER_ONLY)
    get_property(CGAL_LIBRARY TARGET CGAL::CGAL PROPERTY LOCATION)
  else()
    set(CGAL_LIBRARY "")
  endif()
else()
  # We are currently in a CGAL Build and CGALExports.cmake has not
  # necessarily been created yet. Just alias the targets. Also don't
  # access the LOCATION property here to set lib_LIBRARY, since those
  # targets are not imported and this is disallowed by CMP0026. Just
  # set it to the target name.
  if(TARGET CGAL AND NOT TARGET CGAL::CGAL AND NOT CGAL_HEADER_ONLY)
    add_library(CGAL::CGAL ALIAS CGAL)
    set(CGAL_LIBRARY CGAL::CGAL)
  else()
    set(CGAL_LIBRARY "")
  endif()
endif()


# 3RD_PARTY variables.
set(CGAL_3RD_PARTY_INCLUDE_DIRS   "/usr/include")
set(CGAL_3RD_PARTY_DEFINITIONS    "")
set(CGAL_3RD_PARTY_LIBRARIES_DIRS "/usr/lib/x86_64-linux-gnu")
set(CGAL_3RD_PARTY_LIBRARIES      "/usr/lib/x86_64-linux-gnu/libboost_thread.so;/usr/lib/x86_64-linux-gnu/libboost_system.so;/usr/lib/x86_64-linux-gnu/libpthread.so")
