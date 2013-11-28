include(CGAL_suite)
CGAL_example_suite(
  NAME Skin_surface_3_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/NGHK_skin_surface_simple.cpp")
    # The two following .cpp files are in dont_submit, and thus not in the
    # internal releases.
    CGAL_example(SOURCE NGHK_skin_surface_simple.cpp REQUIRED_DEPENDENCIES CGAL)
    CGAL_example(SOURCE NGHK_skin_surface_subdiv.cpp REQUIRED_DEPENDENCIES CGAL)
  endif()

  CGAL_example(SOURCE skin_surface_simple.cpp REQUIRED_DEPENDENCIES CGAL)
  CGAL_example(SOURCE skin_surface_subdiv.cpp REQUIRED_DEPENDENCIES CGAL)
  CGAL_example(SOURCE skin_surface_subdiv_with_normals.cpp REQUIRED_DEPENDENCIES CGAL)
  CGAL_example(SOURCE union_of_balls_simple.cpp REQUIRED_DEPENDENCIES CGAL)
  CGAL_example(SOURCE union_of_balls_subdiv.cpp REQUIRED_DEPENDENCIES CGAL)
  CGAL_example(SOURCE skin_surface_pdb_reader.cpp REQUIRED_DEPENDENCIES CGAL ESBTL)
