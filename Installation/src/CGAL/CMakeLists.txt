include(CGAL_library)

CGAL_define_library(
  NAME CGAL
  BUILD_DEFAULT ON
  REQUIRED_DEPENDENCIES BOOST
  OPTIONAL_DEPENDENCIES GMP MPFR
  SOURCES)

# use_essential_libs()

# include_directories (SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS})

# add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

# link_directories(${CGAL_3RD_PARTY_LIBRARIES_DIRS})

# collect_cgal_library(CGAL "") 

# target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})
