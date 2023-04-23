#ifndef STDAFX_H
#define STDAFX_H

#include <cmath>
#include <cassert>
#include <crtdefs.h>

// STL
#include <algorithm>
#include <map>
#include <vector>
#include <stack>
#include <deque>
#include <fstream>
#include <typeindex>
#include <iterator>

// Windows
#include <windows.h>

// Boost
#include <boost/assert.hpp>
#include <boost/call_traits.hpp>
#include <boost/concept_archetype.hpp>
#include <boost/concept_check.hpp>
#include <boost/config.hpp>
#include <boost/format.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/always.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/apply_wrap.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/clear.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/insert_fwd.hpp>
#include <boost/mpl/iterator_range.hpp>
#include <boost/mpl/iterator_tags.hpp>
#include <boost/mpl/lambda.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/reverse_fold.hpp>
#include <boost/mpl/sequence_tag.hpp>
#include <boost/mpl/set/set0.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/none.hpp>
#include <boost/optional.hpp>
//#include <boost/parameter.hpp>
//#include <boost/parameter/binding.hpp>
//#include <boost/parameter/config.hpp>
//#include <boost/parameter/keyword.hpp>
//#include <boost/parameter/macros.hpp>
//#include <boost/parameter/match.hpp>
//#include <boost/parameter/name.hpp>
//#include <boost/parameter/parameters.hpp>
//#include <boost/parameter/preprocessor.hpp>
//#include <boost/parameter/value_type.hpp>
#include <boost/pending/cstddef.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/comparison/less_equal.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/config/config.hpp>
#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/debug/error.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/intercept.hpp>
#include <boost/preprocessor/facilities/is_empty.hpp>
#include <boost/preprocessor/for.hpp>
#include <boost/preprocessor/identity.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/logical/bool.hpp>
#include <boost/preprocessor/logical/compl.hpp>
#include <boost/preprocessor/logical/not.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repeat.hpp>
#include <boost/preprocessor/repetition/deduce_r.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_shifted.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/selection/max.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/first_n.hpp>
#include <boost/preprocessor/seq/fold_left.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/push_back.hpp>
#include <boost/preprocessor/seq/seq.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/tuple/rem.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <memory>
#include <boost/thread/tss.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/utility.hpp>
#include <boost/variant.hpp>
#include <boost/version.hpp>

// CGAL
//#include <CGAL/AABB_traits.h>
//#include <CGAL/AABB_tree.h>
#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/circulator_bases.h>
//#include <CGAL/Compact_container.h>
#include <CGAL/config.h>
#include <CGAL/Default.h>
#include <CGAL/determinant.h>
#include <CGAL/Dimension.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_2.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_3.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/function_objects.h>
#include <CGAL/Hilbert_policy_tags.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Hilbert_sort_2.h>
#include <CGAL/Hilbert_sort_3.h>
#include <CGAL/Hilbert_sort_base.h>
#include <CGAL/Hilbert_sort_d.h>
#include <CGAL/Hilbert_sort_median_2.h>
#include <CGAL/Hilbert_sort_median_3.h>
#include <CGAL/Hilbert_sort_median_d.h>
#include <CGAL/Hilbert_sort_middle_2.h>
#include <CGAL/Hilbert_sort_middle_3.h>
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/Hilbert_sort_middle_d.h>
#include <CGAL/Image_3.h>
#include <CGAL/TDS_3/internal/Dummy_tds_3.h>
#include <CGAL/Number_types/internal/Exact_type_selector.h>
#include <CGAL/STL_Extension/internal/info_check.h>
//#include <CGAL/internal/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Compare_weighted_squared_radius_3.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Power_side_of_oriented_power_sphere_3.h>
//#include <CGAL/internal/Static_filters/Regular_triangulation_static_filters_traits_3.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
//#include <CGAL/TDS_3/internal/Triangulation_ds_circulators_3.h>
//#include <CGAL/TDS_3/internal/Triangulation_ds_iterators_3.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Kernel/interface_macros.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Lazy.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Lazy_kernel.h>
//#include <CGAL/Mesher_level.h>
//#include <CGAL/Mesher_level_default_implementations.h>
//#include <CGAL/Mesher_level_visitors.h>
//#include <CGAL/Meshes/Filtered_multimap_container.h>
//#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>
//#include <CGAL/Mesh_3/global_parameters.h>
//#include <CGAL/Mesh_3/Mesher_3.h>
//#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>
//#include <CGAL/Mesh_3/mesh_standard_criteria.h>
//#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>
//#include <CGAL/Mesh_3/Mesh_surface_cell_base_3.h>
//#include <CGAL/Mesh_3/Refine_cells_3.h>
//#include <CGAL/Mesh_3/Refine_facets_3.h>
//#include <CGAL/Mesh_3/Refine_tets_visitor.h>
//#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
//#include <CGAL/Mesh_3/Triangle_accessor_primitive.h>
//#include <CGAL/Mesh_3/utilities.h>
//#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>
#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_edge_criteria_3.h>
#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_facet_topology.h>
//#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Multiscale_sort.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Polygon_2/Polygon_2_simplicity.h>
#include <CGAL/Polygon_2/polygon_assertions.h>
#include <CGAL/Polygon_2_algorithms.h>
//#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/predicates/predicates_on_weighted_points_cartesian_3.h>
//#include <CGAL/predicates/Regular_triangulation_ftC3.h>
//#include <CGAL/predicates/Regular_triangulation_rtH3.h>
//#include <CGAL/Profile_counter.h>
//#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Robust_construction.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
//#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/tags.h>
#include <CGAL/Timer.h>
#include <CGAL/Triangle_accessor_3.h>
//#include <CGAL/Triangulation_3.h>
//#include <CGAL/Triangulation_cell_base_3.h>
//#include <CGAL/Triangulation_data_structure_3.h>
//#include <CGAL/Triangulation_ds_cell_base_3.h>
//#include <CGAL/Triangulation_ds_vertex_base_3.h>
//#include <CGAL/Triangulation_simplex_3.h>
//#include <CGAL/Triangulation_structural_filtering_traits.h>
//#include <CGAL/Triangulation_utils_2.h>
//#include <CGAL/Triangulation_utils_3.h>
//#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/tuple.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/utility.h>
//#include <CGAL/Weighted_point.h>

// Mesh_3
/*#include <CGAL_demo/Viewer.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Plugin_helper.h>
#include <ui_Meshing_dialog.h>
#include <Scene_polyhedron_item.h>
#include <implicit_functions/Implicit_function_interface.h>
#include <CGAL_demo/Scene_item_with_display_list.h>*/

#endif //STDAFX_H
