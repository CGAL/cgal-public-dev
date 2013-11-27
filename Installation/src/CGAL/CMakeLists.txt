include(CGAL_library)

CGAL_define_library(
  NAME CGAL
  BUILD_DEFAULT ON
  REQUIRED_DEPENDENCIES Boost Boost_system Boost_thread
  OPTIONAL_DEPENDENCIES GMP MPFR
  SOURCES)
