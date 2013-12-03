if(NOT CGAL_FOUND)
  cmake_minimum_required(VERSION 2.8.11)
  find_package(CGAL REQUIRED)
  include(${CGAL_DIR}/CGAL-export.cmake)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})
endif()

include(CGAL_suite)
CGAL_example_suite(
  NAME Arrangement_on_surface_2_examples
  REQUIRED_DEPENDENCIES CGAL
  OPTIONAL_DEPENDENCIES
  SOURCES aggregated_insertion.cpp algebraic_curves.cpp algebraic_segments.cpp batched_point_location.cpp Bezier_curves.cpp bgl_dual_adapter.cpp bgl_primal_adapter.cpp bounded_planar_vertical_decomposition.cpp circles.cpp circular_arcs.cpp circular_line_arcs.cpp conic_multiplicities.cpp conics.cpp consolidated_curve_data.cpp curve_history.cpp dcel_extension.cpp dcel_extension_io.cpp dual_lines.cpp dual_with_data.cpp edge_insertion.cpp edge_manipulation.cpp edge_manipulation_curve_history.cpp face_extension.cpp face_extension_overlay.cpp generic_curve_data.cpp global_insertion.cpp global_removal.cpp incremental_insertion.cpp io.cpp io_curve_history.cpp io_unbounded.cpp isolated_vertices.cpp observer.cpp overlay.cpp overlay_unbounded.cpp point_location_example.cpp polylines.cpp predefined_kernel.cpp predefined_kernel_non_intersecting.cpp rational_functions.cpp rational_functions_rational_coefficients.cpp special_edge_insertion.cpp spherical_insert.cpp sweep_line.cpp tracing_counting.cpp unbounded_non_intersecting.cpp unbounded_rational_functions.cpp unb_planar_vertical_decomposition.cpp vertical_ray_shooting.cpp)
