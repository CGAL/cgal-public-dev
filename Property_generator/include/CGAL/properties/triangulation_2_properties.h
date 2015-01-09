// #ifndef TRIANGULATION_2_PROPERTIES
// #define TRIANGULATION_2_PROPERTIES

// #include <boost/math/constants/constants.hpp>
// #include <CGAL/properties/tags.h>
// #include <cmath>

// namespace CGAL {
// namespace Properties {
// namespace Triangulation_2 {

// // We give the structure of the documentation grouping here.

// \ingroup PkgProperty_generator_properties

// \defgroup Triangulation_2_properties Triangulation_2
// @{
// \brief A namespace to group together functors computing properties of
// triangulations.

// `#include<CGAL/properties/triangulation_2.h>`

// The functors in this group may be used to compute common
// properties of 2D triangulations. In general, the functions are intended for
// statistics purposes and so we usually return values up to double precision. In
// addition, the functors have been written to take into account the possibility
// that the input may contain infinite edges. Thus by default they will perform
// infinite tests and provide sensible values in these cases. In the case that the
// input is guaranteed not to contain infinite vertices, these tests can be
// disabled by providing the No_finite_test_tag.
// @}


// /// \addtogroup Triangulation_2_properties
// /// @{

// /// \brief Functors taking Triangulation_2 Faces.
// /// \defgroup Face_properties Faces
// /// @{

// /******************************************************************************/

// /**
//  * Functor to compute the area of a Triangulation_2 face.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class area {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef typename Triangulation_2::Face_handle argument_type;
//   /// Type accepted by the functor.
//   typedef double result_type;

//   /// @} \\types
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new area functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   area(const Triangulation_2 &tr) : tr(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the area of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The area of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     if (tr.is_infinite(f))
//       return std::numeric_limits<double>::infinity();

//     // Note, code duplicated below.
//     return CGAL::area(f->vertex(0)->point(), f->vertex(1)->point(),
//                       f->vertex(2)->point());
//   }

//   /// @} //operators
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// struct area<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     // Note, code duplicated above.
//     return CGAL::area(f->vertex(0)->point(), f->vertex(1)->point(),
//                       f->vertex(2)->point());
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute circumradius.
// template <typename Triangulation_2>
// static inline double
// circumradius_internal(const typename Triangulation_2::Face_handle &f) {
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point;
//   const Point &p0 = f->vertex(0)->point();
//   const Point &p1 = f->vertex(1)->point();
//   const Point &p2 = f->vertex(2)->point();
//   // The circumradius is the distance from the centre to any vertex:
//   double r2 = (CGAL::circumcenter(p0, p1, p2) - p0).squared_length();
//   return std::sqrt(r2);
// }
// } // internal namespace

// /**
//  * Functor to compute the circumradius of a Triangulation_2 face.
//  *
//  * We define the circumradius of a face to be the circumradius of the
//  * associated circumcircle of the triangle defined by the face.
//  *
// * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class circumradius {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{
//   typedef typename Triangulation_2::Face_handle argument_type;
//   /// Type returned by the functor.
//   typedef double result_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new circumradius functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   circumradius(const Triangulation_2 &tr) : tr(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the circumradius of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The circumradius of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     // Infinite faces have infinite circumradius.
//     if (tr.is_infinite(f))
//       return std::numeric_limits<double>::infinity();
//     return circumradius_internal<Triangulation_2>(f);
//   }
//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// struct circumradius<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     return circumradius_internal<Triangulation_2>(f);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal function to compute aspect ratio.
// template <typename Triangulation_2>
// static inline double
// aspect_ratio_internal(const typename Triangulation_2::Face_handle &f) {
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;
//   double min = std::numeric_limits<double>::infinity();
//   double max = 0.;
//   for (unsigned short i = 0; i < 3; ++i) {
//     const Point_2 &p1 = f->vertex(i)->point();
//     const Point_2 &p2 = f->vertex((i + 1) % 3)->point();

//     double value = (p1 - p2).squared_length();
//     if (value < min)
//       min = value;
//     if (value > max)
//       max = value;
//   }
//   return std::sqrt(max) / std::sqrt(min);
// }
// }

// /**
//  * Functor to compute the aspect ratio of a Triangulation_2 face.
//  *
//  * We define the aspect ratio of a face to be the ratio of the longest
//  * edge over the shortest edge. We return infinity if the face is infinite.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class aspect_ratio {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef typename Triangulation_2::Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new aspect ratio functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   aspect_ratio(const Triangulation_2 &tr) : tr(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the aspect ratio of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The aspect ratio of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     // Infinite faces have infinite aspect ratio.
//     if (tr.is_infinite(f))
//       return std::numeric_limits<double>::infinity();

//     return aspect_ratio_internal<Triangulation_2>(f);
//   }
//   /// @}
// };

// template <typename Triangulation_2>
// struct aspect_ratio<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     return aspect_ratio_internal<Triangulation_2>(f);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute angle.
// template <typename Triangulation_2>
// static inline double
// angle_internal(const typename Triangulation_2::Face_handle &f,
//                unsigned short i) {
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point;
//   typedef typename Triangulation_2::Geom_traits::Vector_2 Vector;

//   const Point &p1 = f->vertex(i)->point();
//   const Point &p2 = f->vertex((i + 1) % 3)->point();
//   const Point &p3 = f->vertex((i + 2) % 3)->point();

//   Vector v1 = (p2 - p1);
//   Vector v2 = (p3 - p1);

//   v1 = v1 / std::sqrt(v1.squared_length());
//   v2 = v2 / std::sqrt(v2.squared_length());

//   return std::acos((v1 * v2));
// }
// }

// /**
//  * Functor to compute the angle within a  Triangulation_2 face.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class angle {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   typedef typename Triangulation_2::Face_handle Face_handle;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new angle functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   angle(const Triangulation_2 &tr) : tr(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the angle of a face from a Triangulation_2.
//    *
//    * \note This function will return \f$\pi\f$ if the angle is opposite
//    * an infinite vertex, and 0 if the angle is at an infinite vertex.
//    *
//    * \todo Check that this function works correctly.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @param  i The index of the desired angle within the face.
//    * @return   The angle of the given face.
//    */
//   result_type operator()(const Face_handle &f, unsigned short i) const {
//     // There are only three angles in a triangle.
//     CGAL_precondition(0 <= i && i < 3);

//     // Face has an infinite vertex.
//     if (tr.is_infinite(f)) {
//       if (tr.is_infinite(f->vertex(i)))
//         return 0;
//       else
//         return boost::math::constants::pi<double>();
//     }
//     return angle_internal<Triangulation_2>(f, i);
//   }
//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// struct angle<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f, unsigned short i) const {
//     return angle_internal<Triangulation_2>(f, i);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute min_angle.
// template <typename Triangulation_2, typename tag>
// static inline double
// min_angle_internal(const typename Triangulation_2::Face_handle &f,
//                    const angle<Triangulation_2, tag> &angle_functor) {
//   double min = std::numeric_limits<double>::infinity();
//   for (unsigned short i = 0; i < 3; ++i) {
//     double value = angle_functor(f, i);
//     if (value < min)
//       min = value;
//   }

//   return min;
// }
// }

// /**
//  * Functor to compute the minimum angle of a Triangulation_2 face.
//  *
//  * We return the minimum internal angle of the input face, or zero
//  * if one of the vertices is at infinity.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_angle {
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute angles.
//   angle<Triangulation_2> angle_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new minimum angle functor.
//    * @param tr  the Triangulation_ 2 from which the input faces will be
//    *     taken.
//    */
//   min_angle(const Triangulation_2 &tr) : tr(tr), angle_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the minimum angle of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The minimum angle of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     double min = std::numeric_limits<double>::infinity();
//     for (unsigned short i = 0; i < 3; ++i) {
//       double value = angle_functor(f, i);
//       if (value < min)
//         min = value;
//     }
//     return min;
//   }
//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_angle<Triangulation_2, no_finite_test_tag> {
//   // We cache the angle functor locally to save wasted constructions.
//   angle<Triangulation_2, no_finite_test_tag> angle_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     return min_angle_internal<Triangulation_2, no_finite_test_tag>(
//         f, angle_functor);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute max_angle.
// template <typename Triangulation_2, typename tag>
// static inline double
// max_angle_internal(const typename Triangulation_2::Face_handle &f,
//                    const angle<Triangulation_2, tag> &angle_functor) {
//   double max = 0.;
//   for (unsigned short i = 0; i < 3; ++i) {
//     double value = angle_functor(f, i);
//     if (value > max)
//       max = value;
//   }
//   return max;
// }
// }

// /**
//  * Functor to compute the max angle of a Triangulation_2 face.
//  *
//  * Return the maximum internal angle of a face, or \f$\pi\f$ if one
//  * of the vertices is infinite.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_angle {
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute angles.
//   angle<Triangulation_2> angle_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new max angle functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   max_angle(const Triangulation_2 &tr) : tr(tr), angle_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the max angle of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The max angle of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     return max_angle_internal<
//         Triangulation_2, CGAL::Properties::Finite_test_tag>(f, angle_functor);
//   }
//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_angle<Triangulation_2, no_finite_test_tag> {
//   // We cache the angle functor locally to save wasted constructions.
//   angle<Triangulation_2, no_finite_test_tag> angle_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     return max_angle_internal<Triangulation_2, no_finite_test_tag>(
//         f, angle_functor);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute min_edge_length.
// template <typename Triangulation_2, typename tag>
// static inline double
// min_edge_length_internal(const typename Triangulation_2::Face_handle &f,
//                          const Triangulation_2_edge_properties::length<
//                              Triangulation_2, tag> &length_functor) {
//   double min = std::numeric_limits<double>::infinity();

//   for (unsigned short i = 0; i < 3; ++i) {
//     double value = length_functor(f, i);
//     if (value < min)
//       min = value;
//   }
//   return min;
// }
// }

// /**
//  * Functor to compute the min edge length of a Triangulation_2 face.
//  *
//  * We define the minimum edge length as the length of the shortest edge
//  * contained within the input face.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_edge_length {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute edge lengths.
//   Triangulation_2_edge_properties::length<
//       Triangulation_2, CGAL::Properties::Finite_test_tag> length_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef typename Triangulation_2::Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new min edge length functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   min_edge_length(const Triangulation_2 &tr) : tr(tr), length_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the min edge length of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The min edge length of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     return min_edge_length_internal<Triangulation_2, tag>(f, length_functor);
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_edge_length<Triangulation_2, no_finite_test_tag> {
//   Triangulation_2_edge_properties::length<Triangulation_2, no_finite_test_tag>
//   length_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     return min_edge_length_internal<Triangulation_2, no_finite_test_tag>(
//         f, length_functor);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute max_edge_length.
// template <typename Triangulation_2, typename tag>
// static inline double
// max_edge_length_internal(const typename Triangulation_2::Face_handle &f,
//                          const Triangulation_2_edge_properties::length<
//                              Triangulation_2, tag> &length_functor) {
//   double max = 0;
//   for (int i = 0; i < 3; ++i) {
//     double value = length_functor(f, i);
//     if (value < max)
//       max = value;
//   }
//   return max;
// }
// }

// /**
//  * Functor to compute the max edge length of a Triangulation_2 face.
//  *
//  * We define the maximum edge length to be the length of the longest
//  * edge contained within this face, or infinity if one of the vertices
//  * on this face is infinite.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_edge_length {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute edge lengths.
//   Triangulation_2_edge_properties::length<Triangulation_2> length_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef typename Triangulation_2::Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new max edge length functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   max_edge_length(const Triangulation_2 &tr) : tr(tr), length_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the max edge length of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The max edge length of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     return max_edge_length_internal<Triangulation_2, tag>(f, length_functor);
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_edge_length<Triangulation_2, no_finite_test_tag> {
//   Triangulation_2_edge_properties::length<Triangulation_2, no_finite_test_tag>
//   length_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     // return aspect_ratio_internal(f);
//     return max_edge_length_internal<Triangulation_2, no_finite_test_tag>(
//         f, length_functor);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute min_neighbor_area.
// template <typename Triangulation_2, typename tag>
// static inline double
// min_neighbor_area_internal(typename Triangulation_2::Face_handle f,
//                            const area<Triangulation_2, tag> &area_functor) {
//   double min = std::numeric_limits<double>::infinity();
//   for (unsigned short i = 0; i < 3; ++i) {
//     double value = area_functor(f->neighbor(i));
//     if (value < min)
//       min = value;
//   }
//   return min;
// }
// }

// /**
//  * Functor to compute the min neighbor area of a Triangulation_2 face.
//  *
//  * We define the minimum neighbor area to be the smallest area of any
//  * face sharing an edge with this face.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_neighbor_area {
//   typedef typename Triangulation_2::Face_handle Face_handle;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute angles.
//   area<Triangulation_2> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new min neighbor area functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   min_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the min neighbor area of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The min neighbor area of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     return min_neighbor_area_internal<
//         Triangulation_2, CGAL::Properties::Finite_test_tag>(f, area_functor);
//     ;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   area<Triangulation_2, no_finite_test_tag> area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     // return aspect_ratio_internal(f);
//     return min_neighbor_area_internal<Triangulation_2, no_finite_test_tag>(
//         f, area_functor);
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute max_neighbor_area_internal.
// template <typename Triangulation_2, typename tag>
// static inline double
// max_neighbor_area_internal(typename Triangulation_2::Face_handle f,
//                            const area<Triangulation_2, tag> &area_functor) {
//   double max = 0;
//   for (unsigned short i = 0; i < 3; ++i) {
//     double value = area_functor(f->neighbor(i));
//     if (value > max)
//       max = value;
//   }
//   return max;
// }
// }

// /**
//  * Functor to compute the max neighbor area of a Triangulation_2 face.
//  *
//  * We define the maximum neighbor area to be the area of the face
//  * with the greatest area sharing an edge with this face, or infinity
//  * if one of them contains the infinite vertex.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_neighbor_area {
//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute angles.
//   area<Triangulation_2> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;
//   /// Type accepted by the functor.
//   typedef typename Triangulation_2::Face_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new max neighbor area functor.
//    * @param tr  the Triangulation_2 from which the input faces will be
//    *     taken.
//    */
//   max_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute the max neighbor area of a face from a Triangulation_2.
//    *
//    * \pre
//    * The face is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f The face.
//    * @return   The max neighbor area of the given face.
//    */
//   result_type operator()(const argument_type &f) const {
//     return max_neighbor_area_internal<
//         Triangulation_2, CGAL::Properties::Finite_test_tag>(f, area_functor);
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   area<Triangulation_2, no_finite_test_tag> area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Face_handle argument_type;
//   double operator()(const argument_type &f) const {
//     // return aspect_ratio_internal(f);
//     return max_neighbor_area_internal<Triangulation_2, no_finite_test_tag>(
//         f, area_functor);
//   }
// };

// /******************************************************************************/
// // Helper functions.
// /******************************************************************************/

// /**
//  * Construct a functor to compute the XYZ for Triangulation_2 Face
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// area<Triangulation_2, tag> make_area_functor(const Triangulation_2 &tr_2,
//                                               tag) {
//   return area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the circumradius for Triangulation_2 Face
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// circumradius<Triangulation_2, tag>
// make_circumradius_functor(const Triangulation_2 &tr_2, tag) {
//   return circumradius<Triangulation_2>(tr_2);
// }

// /**
//  * Construct a functor to compute the aspect_ratio for Triangulation_2 Face
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// aspect_ratio<Triangulation_2, tag>
// make_aspect_ratio_functor(const Triangulation_2 &tr_2, tag) {
//   return aspect_ratio<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the min_edge_length for Triangulation_2
//  * Face handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// min_edge_length<Triangulation_2, tag>
// make_min_edge_length_functor(const Triangulation_2 &tr_2, tag) {
//   return min_edge_length<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the max_edge_length for Triangulation_2
//  * Face handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_edge_length<Triangulation_2, tag>
// make_max_edge_length_functor(const Triangulation_2 &tr_2, tag) {
//   return max_edge_length<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the min_angle for Triangulation_2 Face
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// min_angle<Triangulation_2, tag>
// make_min_angle_functor(const Triangulation_2 &tr_2, tag) {
//   return min_angle<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the max_angle for Triangulation_2 Face
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_angle<Triangulation_2, tag>
// make_max_angle_functor(const Triangulation_2 &tr_2, tag) {
//   return max_angle<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the min_neighbor_area for Triangulation_2
//  * Face handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// min_neighbor_area<Triangulation_2, tag>
// make_min_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return min_neighbor_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the max_neighbor_area for Triangulation_2
//  * Face handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_neighbor_area<Triangulation_2, tag>
// make_max_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return max_neighbor_area<Triangulation_2>(tr_2, tag);
// }

// /// @} // End of Edge_properties

// /*******************************************************************************
// *
// *
// *
// *
// *   We now give properties for triangulation vertices.
// *
// *
// *
// *
// *******************************************************************************/

// /// \brief Functors taking Triangulation_2 Vertices.
// /// \defgroup Vertex_properties Vertices
// /// @{

// /******************************************************************************/

// /**
//  * Functor to compute the degree of a `Triangulation_2` `Vertex_handle`.
//  * This provides a wrapper to the function provided within the
//  * vertex handle class for consistancy.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2> class degree {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef int result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new degree functor.
//    */
//   degree() {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute degree of a `Vertex_handle` from a `Triangulation_2.`
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The degree of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const { return v->degree(); }

//   /// @}
// };

// /******************************************************************************/

// /**
//  * Functor to compute the area of the Voronoi cell associated with this
//  * vetrex.
//  *
//  * \todo Check that this definition will work with any Triangulation_2 type.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class dual_area {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Geom_traits Gt;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;

//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new dual area functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   dual_area(const Triangulation_2 &tr) : tr(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute dual area of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The dual area of the given vertex.
//    */
//   // NOTE: Code duplicated below.
//   result_type operator()(const argument_type &v) const {
//     Polygon_2<Gt> p;

//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       if (tr.is_infinite(f))
//         return std::numeric_limits<double>::infinity();
//       do
//         p.push_back(tr.circumcenter(f));
//       while (++f != done);
//     }

//     return fabs(p.area());
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class dual_area<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;

//   // We need a reference to the triangulation to compute this property.
//   dual_area(const Triangulation_2 &tr) : tr(tr) {}

//   // NOTE: Code duplicated above.
//   double operator()(const argument_type &v) const {
//     Polygon_2<typename Triangulation_2::Geom_traits> p;
//     Face_circulator c = tr.incident_faces(v), done(c);
//     if (c != 0) {
//       do
//         p.push_back(tr.circumcenter(c));
//       while (++c != done);
//     }

//     return fabs(p.area());
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the star area of this vetrex.
//  *
//  * We define the star area to be the sum of the areas of the triangles
//  * that contain this vertex, or infinity if one of the vertices is
//  * infinite.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class star_area {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;

//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::area<Triangulation_2, tag> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new star area functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   star_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute star area of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The star area of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const {
//     // NOTE: Code duplicated below.
//     double area = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         area += area_functor(f);
//       } while (++f != done);
//     }
//     return area;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class star_area<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::area<Triangulation_2, no_finite_test_tag>
//   area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   star_area(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double area = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         area += area_functor(f);
//       } while (++f != done);
//     }
//     return area;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the link length of this vetrex.
//  *
//  * We define the link length to be the sum of the lengths of all edges not
//  * containing this point, but contained within all triangles containing this
//  * vertex.
//  *
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class link_length {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;

//   const Triangulation_2 &tr;
//   Triangulation_2_edge_properties::length<Triangulation_2> length_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new link length functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   link_length(const Triangulation_2 &tr) : tr(tr), length_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute link length of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The link length of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const {
//     // NOTE: Code duplicated below.
//     double length = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do
//         length += length_functor(f, f->index(v));
//       while (++f != done);
//     }
//     return length;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class link_length<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_edge_properties::length<Triangulation_2, no_finite_test_tag>
//   length_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   link_length(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double length = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do
//         length += length_functor(f, f->index(v));
//       while (++f != done);
//     }
//     return length;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the max neighbour area of this vetrex.
//  *
//  * We define the maximum neighbor area to be the area of the largest
//  * face neighboring this vertex, or infinity if one of the faces is infinte.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_neighbor_area {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;

//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::area<Triangulation_2> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new max neighbour area functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   max_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute max neighbour area of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The max neighbour area of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const {
//     // Note: Code duplicated below.
//     double max = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = area_functor(f);
//         if (value > max)
//           max = value;
//       } while (++f != done);
//     }
//     return max;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::area<Triangulation_2, no_finite_test_tag>
//   area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   max_neighbor_area(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double max = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = area_functor(f);
//         if (value > max)
//           max = value;
//       } while (++f != done);
//     }
//     return max;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the minimum neighbor area of this vetrex.
//  *
//  * We define the minum neighbour area to be the area of the smallest
//  * triangle neighboring the input vertex.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_neighbor_area {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;
//   const Triangulation_2 &tr;
//   typename Triangulation_2_face_properties::area<Triangulation_2> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new minimum neighbor area functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *            taken.
//    */
//   min_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute minimum neighbor area of a Vertex_handle from a
//    * Triangulation_2.
//    *

//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The minimum neighbor area of the given vertex.
//    */

//   result_type operator()(const argument_type &v) const {
//     // NOTE: Code duplicated below.
//     double min = std::numeric_limits<double>::infinity();
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = area_functor(f);
//         if (value < min)
//           min = value;
//       } while (++f != done);
//     }
//     return min;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::area<Triangulation_2, no_finite_test_tag>
//   area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   min_neighbor_area(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double min = std::numeric_limits<double>::infinity();

//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = area_functor(f);
//         if (value < min)
//           min = value;
//       } while (++f != done);
//     }
//     return min;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the max star angle of this vetrex.
//  *
//  * We define the maximum star angle to be the largest angle defined by any
//  * two edges exiting this vertex.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_star_angle {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::angle<Triangulation_2> angle_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new max star angle functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   max_star_angle(const Triangulation_2 &tr) : tr(tr), angle_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute max star angle of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The max star angle of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const {
//     // NOTE: Code duplicated below.
//     double max = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = angle_functor(f, f->index(v));
//         if (value > max)
//           max = value;
//       } while (++f != done);
//     }
//     return max;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_star_angle<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::angle<Triangulation_2, no_finite_test_tag>
//   angle_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   max_star_angle(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double max = 0.;
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = angle_functor(f, f->index(v));
//         if (value > max)
//           max = value;
//       } while (++f != done);
//     }
//     return max;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the min star angle of this vetrex.
//  *
//  * We define the min star angle to be the smallest angle defined between
//  * any two edges exiting this vertex.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_star_angle {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_circulator Face_circulator;
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::angle<Triangulation_2, tag> angle_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Vertex_handle argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new min star angle functor.
//    * @param tr  the Triangulation_2 from which the Vertex handles will be
//    *     taken.
//    */
//   min_star_angle(const Triangulation_2 &tr) : tr(tr), angle_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute min star angle of a Vertex_handle from a Triangulation_2.
//    *
//    * \pre
//    * The vertex is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  v The input Vertex_handle.
//    * @return   The min star angle of the given vertex.
//    */
//   result_type operator()(const argument_type &v) const {
//     // NOTE: Code duplicated below.
//     double min = std::numeric_limits<double>::infinity();
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = angle_functor(f, f->index(v));
//         if (value < min)
//           min = value;
//       } while (++f != done);
//     }
//     return min;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_star_angle<Triangulation_2, no_finite_test_tag> {
//   const Triangulation_2 &tr;
//   Triangulation_2_face_properties::angle<Triangulation_2, no_finite_test_tag>
//   angle_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Vertex_handle argument_type;
//   min_star_angle(const Triangulation_2 &tr) : tr(tr) {}
//   double operator()(const argument_type &v) const {
//     // NOTE: Code duplicated above.
//     double min = std::numeric_limits<double>::infinity();
//     Face_circulator f = tr.incident_faces(v), done(f);
//     if (f != 0) {
//       do {
//         double value = angle_functor(f, f->index(v));
//         if (value < min)
//           min = value;
//       } while (++f != done);
//     }
//     return min;
//   }
// };

// /******************************************************************************/

// /**
//  * Construct a functor to compute the min_star_angle for Triangulation_2
//  * Vertex handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// min_star_angle<Triangulation_2, tag>
// make_min_star_angle_functor(const Triangulation_2 &tr_2) {
//   return min_star_angle<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the max_star_angle for Triangulation_2
//  * Vertex handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_star_angle<Triangulation_2, tag>
// make_max_star_angle_functor(const Triangulation_2 &tr_2, tag) {
//   return max_star_angle<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the min_neighbor_area for Triangulation_2
//  * Vertex handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// min_neighbor_area<Triangulation_2, tag>
// make_min_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return min_neighbor_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the max_neighbor_area for Triangulation_2
//  * Vertex handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_neighbor_area<Triangulation_2, tag>
// make_max_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return max_neighbor_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the link_length for Triangulation_2 Vertex
//  * handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// link_length<Triangulation_2, tag>
// make_link_length_functor(const Triangulation_2 &tr_2, tag) {
//   return link_length<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the star_area for Triangulation_2 Vertex
//  * handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// star_area<Triangulation_2, tag>
// make_star_area_functor(const Triangulation_2 &tr_2, tag) {
//   return star_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the dual_area for Triangulation_2 Vertex
//  * handles.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// dual_area<Triangulation_2, tag>
// make_dual_area_functor(const Triangulation_2 &tr_2, tag) {
//   return dual_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the degree for Triangulation_2 Vertex
//  * handles.
//  *
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// degree<Triangulation_2, tag> make_degree_functor(const Triangulation_2 &tr_2,
//                                                   tag) {
//   return degree<Triangulation_2, tag>();
// }

// /// @}

// /*******************************************************************************
// *
// *
// *
// *
// *   We now give properties for triangulation Faces.
// *
// *
// *
// *
// *******************************************************************************/

// /// \brief Functors taking Triangulation_2 Edges.
// /// \defgroup Edge_properties Edges
// /// @{

// /******************************************************************************/

// /**
//  * Functor to compute the length of a Triangulation_2 edge.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class length {
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Edge Edge;
//   typedef typename Triangulation_2::Face_handle Face_handle;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Edge argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new length functor.
//    * @param tr  the Triangulation_2 from which the edges will be taken.
//    */
//   length(const Triangulation_2 &tr) : tr(tr) {};

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute length of an Edge from a Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  e The edge.
//    * @return   The length of the given edge or infinity if the edge is not
//    *           finite.
//    */
//   result_type operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }

//   /**
//    * Functor operator to compute length of an Edge from a Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f Face_handle containing our desired edge.
//    * @param  i Index of the edge in the given face.
//    * @return   The length of the given edge or infinity if the edge is not
//    *           finite.
//    */
//   result_type operator()(const Face_handle &f, unsigned short i) const {
//     const Vertex_handle &v1 = f->vertex(f->cw(i));
//     const Vertex_handle &v2 = f->vertex(f->ccw(i));

//     if (tr.is_infinite(v1))
//       return std::numeric_limits<double>::infinity();
//     if (tr.is_infinite(v2))
//       return std::numeric_limits<double>::infinity();

//     // Don't use an internal function for something this simple.
//     // We just duplicate the code.
//     return std::sqrt((v1->point() - v2->point()).squared_length());
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// struct length<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Edge argument_type;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   double operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }
//   double operator()(const Face_handle &f, unsigned short i) const {
//     // Don't use an internal function for something this simple.
//     // We just duplicate the code.
//     const Vertex_handle &v1 = f->vertex(f->cw(i));
//     const Vertex_handle &v2 = f->vertex(f->ccw(i));
//     return std::sqrt((v1->point() - v2->point()).squared_length());
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the neighbor area of a Triangulation_2 Edge.
//  *
//  * We define the 'neighbor area' to be the sum of the area of the two faces
//  * adjoining a given edge.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class neighbor_area {
//   typedef typename Triangulation_2::Edge Edge;
//   typedef typename Triangulation_2::Face_handle Face_handle;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute face area.
//   Triangulation_2_face_properties::area<
//       Triangulation_2, CGAL::Properties::Finite_test_tag> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Edge argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new neighbor area functor.
//    * @param tr  the Triangulation_2 from which the edges will be taken.
//    */
//   neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute neighbor area of an Edge from a Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  e The edge.
//    * @return   The neighbor area of the given edge or infinity if one of
//    *           the adjoining faces is not finite.
//    */
//   result_type operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }

//   /**
//    * Functor operator to compute the neighbor area of an Edge from a
//    * Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f Face_handle containing our desired edge.
//    * @param  i Index of the edge in the given face.
//    * @return   The neighbor area of the given edge or infinity if one of
//    *           the neighboring points is not finite.
//    */
//   result_type operator()(const Face_handle &f, int i) const {
//     // NOTE: code is duplicated below.
//     return area_functor(f) + area_functor(f->neighbor(i));
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class neighbor_area<Triangulation_2, no_finite_test_tag> {
//   Triangulation_2_face_properties::area<Triangulation_2, no_finite_test_tag>
//   area_functor;
//   typedef typename Triangulation_2::Face_handle Face_handle;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Edge argument_type;
//   double operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }
//   double operator()(const Face_handle &f, unsigned short i) const {
//     // NOTE: code is duplicated above.
//     return area_functor(f) + area_functor(f->neighbor(i));
//   }
// };

// /******************************************************************************/

// // Limit this function to use within this unit alone.
// namespace {
// // Internal funtion to compute circumradius.
// // We take two faces, as it makes the code slightly more efficient.
// template <typename Triangulation_2>
// static inline double
// dual_length_internal(const typename Triangulation_2::Face_handle &f1,
//                      const typename Triangulation_2::Face_handle &f2) {
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point;

//   const Point &p1 = f1->vertex(0)->point();
//   const Point &p2 = f1->vertex(1)->point();
//   const Point &p3 = f1->vertex(2)->point();

//   const Point &q1 = f2->vertex(0)->point();
//   const Point &q2 = f2->vertex(1)->point();
//   const Point &q3 = f2->vertex(2)->point();

//   const Point &c1 = CGAL::circumcenter(p1, p2, p3);
//   const Point &c2 = CGAL::circumcenter(q1, q2, q3);
//   return std::sqrt((c1 - c2).squared_length());
// }
// }

// /**
//  * Functor to compute the dual length of a Triangulation_2 Edge.
//  *
//  * We define the 'dual length' to be the length of the Voronoi edge
//  * between the two adjacent triangles.
//  *
//  * \todo Check the dual edge definition is sensible for all
//  *        possible input classes.
//  *
//  * \note
//  * Iterating over finite edges in a Triangulation_2 will always result in
//  * infinitely long dual edges.
//  *
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class dual_length {
//   typedef typename Triangulation_2::Edge Edge;
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Edge argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new neighbor area functor.
//    * @param tr  the Triangulation_2 from which the edges will be taken.
//    */
//   dual_length(const Triangulation_2 &tr) : tr(tr) {};

//   /// @}
//   /// \name Operators
//   /// @{

//   **Compute dual length of an Edge from a Triangulation_2. ** \pre *The edge is
//       contained within the Triangulation_2 provided on *construction. **@param e
//       The edge. *@ return The dual length of the given edge or infinity if one *
//       of the dual vertices is
//       infinity.result_type operator()(const argument_type & e) const {
//     const Face_handle &f1 = e.first;
//     const Face_handle &f2 = f1->neighbor(e.second);

//     if (tr.is_infinite(f1) || tr.is_infinite(f2))
//       return std::numeric_limits<double>::infinity();

//     return dual_length_internal<Triangulation_2>(f1, f2);
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// struct dual_length<Triangulation_2, no_finite_test_tag> {
//   typedef double result_type;
//   typedef typename Triangulation_2::Edge argument_type;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   double operator()(const argument_type &e) const {
//     const Face_handle &f1 = e.first;
//     const Face_handle &f2 = f1->neighbor(e.second);
//     return dual_length_internal<Triangulation_2>(f1, f2);
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the max neighbor area of a Triangulation_2 Edge.
//  *
//  * We define the 'max neighbor area' to be the area of the adjacent triangle
//  * with the largest area, or infinity if one of the triangles is infinite.
//  *
//  * \cgalAdvancedBegin
//  *   This functor maintains a reference to the Triangulation_2
//  *   it is given on construction. This allows the functor to check
//  *   whether or not a given input contains any points which are infinite.
//  * \cgalAdvancedEnd
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class max_neighbor_area {
//   typedef typename Triangulation_2::Edge Edge;
//   typedef typename Triangulation_2::Vertex_handle Vertex_handle;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute face area.
//   typename Triangulation_2_face_properties::area<
//       Triangulation_2, CGAL::Properties::Finite_test_tag> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Edge argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new neighbor area functor.
//    * @param tr  the Triangulation_2 from which the edges will be taken.
//    */
//   max_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute max neighbor area of an Edge from a Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  e The edge.
//    * @return   The max neighbor area of the given edge.
//    */
//   result_type operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }

//   // Given a Face_handle and an index.
//   result_type operator()(const Face_handle &f, unsigned short i) const {
//     // NOTE: Code duplicated below.
//     const Face_handle &f2 = f->neighbor(i);
//     double area_1 = area_functor(f);
//     double area_2 = area_functor(f2);
//     if (area_1 > area_2)
//       return area_1;
//     return area_2;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class max_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   typename Triangulation_2_face_properties::area<
//       Triangulation_2, no_finite_test_tag> area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Edge argument_type;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   double operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }
//   double operator()(const Face_handle &f, unsigned short i) const {
//     // NOTE: Code duplicated above.
//     const Face_handle &f2 = f->neighbor(i);
//     double area_1 = area_functor(f);
//     double area_2 = area_functor(f2);
//     if (area_1 > area_2)
//       return area_1;
//     return area_2;
//   }
// };

// /******************************************************************************/

// /**
//  * Functor to compute the min neighbor area of a Triangulation_2 Edge.
//  *
//  * We define the 'min neighbor area' to be the area of the adjacent triangle
//  * having the smallest area.
//  *
//  * By default the functor checks for finite vertices, but this checking
//  * can be removed by supplying the tag `CGAL::No_finite_test` tag.
//  *
//  * @tparam Triangulation_2 The Triangulation type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// class min_neighbor_area {

//   typedef typename Triangulation_2::Edge Edge;
//   typedef typename Triangulation_2::Face_handle Face_handle;

//   // Allow us to keep a reference to the triangulation.
//   const Triangulation_2 &tr;

//   // Functor to compute face area.
//   typename Triangulation_2_face_properties::area<
//       Triangulation_2, CGAL::Properties::Finite_test_tag> area_functor;

// public:
//   /// \name Types
//   /// @{

//   /// Type returned by the functor.
//   typedef double result_type;

//   /// Type accepted by the functor.
//   typedef Edge argument_type;

//   /// @}
//   /// \name Constructor
//   /// @{

//   /**
//    * Construct a new neighbor area functor.
//    * @param tr  the Triangulation_2 from which the edges will be taken.
//    */
//   min_neighbor_area(const Triangulation_2 &tr) : tr(tr), area_functor(tr) {}

//   /// @}
//   /// \name Operators
//   /// @{

//   /**
//    * Compute min neighbor area of an Edge from a Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  e The edge.
//    * @return   The min neighbor area of the given edge.
//    */
//   result_type operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }

//   /**
//    * Operator to compute the min neighbor area of an Edge from a
//    * Triangulation_2.
//    *
//    * \pre
//    * The edge is contained within the Triangulation_2 provided on
//    * construction.
//    *
//    * @param  f Face_handle containing our desired edge.
//    * @param  i  Index of the edge in the given face.
//    * @return    The min neighbor area of the given edge.
//    */
//   result_type operator()(const Face_handle &f, unsigned short i) const {
//     // NOTE: Code duplicated below.
//     const Face_handle &f2 = f->neighbor(i);
//     double area_1 = area_functor(f);
//     double area_2 = area_functor(f2);
//     if (area_1 < area_2)
//       return area_1;
//     return area_2;
//   }

//   /// @}
// };

// // Specialisation to disable finiteness tests.
// template <typename Triangulation_2>
// class min_neighbor_area<Triangulation_2, no_finite_test_tag> {
//   typename Triangulation_2_face_properties::area<
//       Triangulation_2, no_finite_test_tag> area_functor;

// public:
//   typedef double result_type;
//   typedef typename Triangulation_2::Edge argument_type;
//   typedef typename Triangulation_2::Face_handle Face_handle;
//   double operator()(const argument_type &e) const {
//     return operator()(e.first, e.second);
//   }
//   double operator()(const Face_handle &f, unsigned short i) const {
//     // NOTE: Code duplicated above.
//     const Face_handle &f2 = f->neighbor(i);
//     double area_1 = area_functor(f);
//     double area_2 = area_functor(f2);
//     if (area_1 < area_2)
//       return area_1;
//     return area_2;
//   }
// };

// /******************************************************************************/
// /******************************************************************************/

// /**
//  * Construct a functor to compute edge lengths.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// length<Triangulation_2, tag> make_length_functor(const Triangulation_2 &tr_2,
//                                                   tag) {
//   return length<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the neighbor area of an edge.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// neighbor_area<Triangulation_2, tag>
// make_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return neighbor_area<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the dual length of an edge.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// dual_length<Triangulation_2, tag>
// make_dual_length_functor(const Triangulation_2 &tr_2, tag) {
//   return dual_length<Triangulation_2, tag>(tr_2);
// }

// /**
//  * Construct a functor to compute the maximum neighbor area for edges.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2,
//           typename tag = CGAL::Properties::Finite_test_tag>
// max_neighbor_area<Triangulation_2, tag>
// make_max_neighbor_area_functor(const Triangulation_2 &tr_2, tag) {
//   return max_neighbor_area<Triangulation_2>(tr_2, tag);
// }

// /**
//  * Construct a functor to compute minimum neighbor areas for edges.
//  *
//  * @param  tr_2            The Triangulation_2 to be associated with the
//  *                         functor.
//  * @tparam Triangulation_2 The Triangulation_2 type.
//  */
// template <typename Triangulation_2>
// min_neighbor_area<Triangulation_2, tag>
// make_min_neighbor_area_functor(const Triangulation_2 &tr_2) {
//   return min_neighbor_area<Triangulation_2>(tr_2, tag);
// }

// /// @} // Vertex properties.

// /******************************************************************************/

// /// @}

// } // namespace Triangulation_2
// } // namespace Properties
// } // namespace CGAL

// /******************************************************************************/
// #endif
// /******************************************************************************/
