// TODO: Copyright info

#ifndef CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
#define CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/aff_transformation_tags.h>

#include <boost/type_traits/is_same.hpp>

#include <pointmatcher/PointMatcher.h>

#include <iostream>
#include <string>
#include <map>

namespace CGAL {

namespace pointmatcher {

template<typename Scalar>
using ICP = typename PointMatcher<Scalar>::ICP;

struct ICP_config {
  std::string name;
  std::map<std::string, std::string> params;
};

namespace internal {

void dump_invalid_point_matcher_config_exception_msg(const PointMatcherSupport::InvalidElement& err) {
  std::cerr << "ERROR Invalid configuration for PM::ICP, omitting configuration: " << std::endl;
	std::cerr << "   " << err.what() << std::endl;
}

template<typename Scalar, typename NamedParameters1, typename NamedParameters2>
ICP<Scalar>
construct_icp(const NamedParameters1& np1, const NamedParameters2& np2)
{
  typedef PointMatcher<Scalar> PM;

  using boost::choose_param;
  using boost::get_param;
  
  ICP<Scalar> icp;

  icp.setDefault();

  ICP_config null_config { .name = "_null_pm_config_in_cgal" };
  auto is_null_config = [&](const ICP_config& c) { return !c.name.compare(null_config.name); };

  // In CGAL, point_set_1 is the reference while point_set_2 is the data
  // However, in pointmatcher, the order is reverse: point_set_1 is the data while point_set_2 is the reference
  // Therefore, filter params from np1 applies to reference data points while params from np2 applies to reading data points

  // np1.point_set_filters -> PM::ReferenceDataPointsFilter
  auto reference_data_points_filter_configs = choose_param(get_param(np1, internal_np::point_set_filters), std::vector<ICP_config>());
  for(const auto& conf : reference_data_points_filter_configs)
  {
    try {
      icp.referenceDataPointsFilters.push_back( PM::get().DataPointsFilterRegistrar.create(conf.name, conf.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // np2.point_set_filters -> PM::ReadingDataPointsFilters
  auto reading_data_points_filter_configs = choose_param(get_param(np2, internal_np::point_set_filters), std::vector<ICP_config>());
  for(const auto& conf : reading_data_points_filter_configs)
  {
    try {
      icp.readingDataPointsFilters.push_back( PM::get().DataPointsFilterRegistrar.create(conf.name, conf.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Matcher
  auto matcher_config = choose_param(get_param(np1, internal_np::matcher), null_config);
  if(!is_null_config(matcher_config))
  {
    try {
      icp.matcher = PM::get().MatcherRegistrar.create(matcher_config.name, matcher_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Outlier Filters
  auto outlier_filters_config = choose_param(get_param(np1, internal_np::outlier_filters), std::vector<ICP_config>());
  for(const auto& conf : outlier_filters_config)
  {
    try {
      icp.outlierFilters.push_back( PM::get().OutlierFilterRegistrar.create(conf.name, conf.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Error Minimizer
  auto error_minimizer_config = choose_param(get_param(np1, internal_np::error_minimizer), null_config);
  if(!is_null_config(error_minimizer_config))
  {
    try {
      icp.errorMinimizer = PM::get().ErrorMinimizerRegistrar.create(error_minimizer_config.name, error_minimizer_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Transformation Checkers
  auto transformation_checkers_config = choose_param(get_param(np1, internal_np::transformation_checkers), std::vector<ICP_config>());
  for(const auto& conf : transformation_checkers_config)
  {
    try {
      icp.transformationCheckers.push_back( PM::get().TransformationCheckerRegistrar.create(conf.name, conf.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Inspector
  auto inspector_config = choose_param(get_param(np1, internal_np::inspector), null_config);
  if(!is_null_config(error_minimizer_config))
  {
    try {
      icp.inspector = PM::get().InspectorRegistrar.create(inspector_config.name, inspector_config.params);
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  // Logger
  auto logger_config = choose_param(get_param(np1, internal_np::logger), null_config);
  if(!is_null_config(logger_config))
  {
    try {
      PointMatcherSupport::setLogger( PM::get().LoggerRegistrar.create(logger_config.name, logger_config.params) );
    } catch(typename PointMatcherSupport::InvalidElement& error) {
      dump_invalid_point_matcher_config_exception_msg(error);
    }
  }

  return icp;
}

template<typename Scalar,
         typename PointRange,
         typename PointMap,
         typename VectorMap,
         typename PM_matrix>
void
copy_cgal_points_to_pm_matrix
(const PointRange& prange, PointMap point_map, VectorMap vector_map, PM_matrix& pm_points, PM_matrix& pm_normals)
{
  int idx = 0;
  for(const auto& p : prange)
  {
    // position
    const auto& pos = get(point_map, p);
    pm_points(0, idx) = pos.x();
    pm_points(1, idx) = pos.y();
    pm_points(2, idx) = pos.z();
    pm_points(3, idx) = Scalar(1.);
    
    // normal
    const auto& normal = get (vector_map, p);
    pm_normals(0, idx) = normal.x();
    pm_normals(1, idx) = normal.y();
    pm_normals(2, idx) = normal.z();
    
    ++idx;
  }
}

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
typename Kernel::Aff_transformation_3
compute_registration_transformation(const PointRange1& range1, const PointRange2& range2,
                                    PointMap1 point_map1, PointMap2 point_map2,
                                    VectorMap1 vector_map1, VectorMap2 vector_map2,
                                    const typename Kernel::Aff_transformation_3& initial_transform,
                                    ICP<typename Kernel::FT> icp)
{
  using Scalar    = typename Kernel::FT;
  
  using PM                  = PointMatcher<Scalar>;
  using PM_cloud            = typename PM::DataPoints;
  using PM_matrix           = typename PM::Matrix;
  using PM_labels           = typename PM_cloud::Labels;
  using PM_transform        = typename PM::Transformation;
  using PM_transform_params = typename PM::TransformationParameters;
  
  // ref_points: 1, points: 2
  std::size_t nb_ref_points = range1.size();
  std::size_t nb_points     = range2.size();
  
  PM_matrix ref_points_pos_matrix    = PM_matrix (4, nb_ref_points);
  PM_matrix ref_points_normal_matrix = PM_matrix (3, nb_ref_points);
  PM_matrix points_pos_matrix    = PM_matrix (4, nb_points);
  PM_matrix points_normal_matrix = PM_matrix (3, nb_points);

  // In CGAL, point_set_1 is the reference while point_set_2 is the data

  // convert cgal points to pointmatcher points
  internal::copy_cgal_points_to_pm_matrix<Scalar>(range1,
                                                  point_map1,
                                                  vector_map1,
                                                  ref_points_pos_matrix, // out
                                                  ref_points_normal_matrix); // out

  internal::copy_cgal_points_to_pm_matrix<Scalar>(range2,
                                                  point_map2,
                                                  vector_map2,
                                                  points_pos_matrix, // out
                                                  points_normal_matrix); // out
                                        
  auto construct_PM_cloud = [](const PM_matrix& positions, const PM_matrix& normals) -> PM_cloud
  {
    PM_cloud cloud;
    
    cloud.addFeature("x", positions.row(0));
    cloud.addFeature("y", positions.row(1));
    cloud.addFeature("z", positions.row(2));
    cloud.addFeature("pad", positions.row(3));
    cloud.addDescriptor("normals", normals);
    
    return cloud;
  };

  PM_cloud ref_cloud = construct_PM_cloud(ref_points_pos_matrix, ref_points_normal_matrix);
  PM_cloud cloud     = construct_PM_cloud(points_pos_matrix,     points_normal_matrix);
  
  PM_transform_params pm_transform_params = PM_transform_params::Identity(4,4);

  // Convert CGAL transform to pm transform
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      pm_transform_params(i,j) = initial_transform.m(i,j);

  try 
  {
		const PM_transform_params prior = pm_transform_params;
		pm_transform_params = icp(cloud, ref_cloud, prior);
	}
	catch (typename PM::ConvergenceError& error) // TODO: Shall we make it a CGAL exception?
	{
		std::cerr << "ERROR PM::ICP failed to converge: " << std::endl;
		std::cerr << "   " << error.what() << std::endl;
	}

	// Rigid transformation
	std::shared_ptr<PM_transform> transform = PM::get().REG(Transformation).create("RigidTransformation");
  pm_transform_params = transform->correctParameters(pm_transform_params);

  typename Kernel::Aff_transformation_3 cgal_transform
    (pm_transform_params(0,0), pm_transform_params(0,1), pm_transform_params(0,2), pm_transform_params(0,3),
     pm_transform_params(1,0), pm_transform_params(1,1), pm_transform_params(1,2), pm_transform_params(1,3),
     pm_transform_params(2,0), pm_transform_params(2,1), pm_transform_params(2,2), pm_transform_params(2,3));

#ifdef CGAL_POINTMATCHER_VERBOSE
  std::cerr << "Transformation matrix: " << std::endl;
  for (std::size_t i = 0; i < 4; ++ i)
  {
    for (std::size_t j = 0; j < 4; ++ j)
      std::cerr << cgal_transform.coeff(i,j) << " ";
    std::cerr << std::endl;
  }
#endif

  return cgal_transform;
}

} // end of namespace internal

// TODO: Document
// point_set_1 is reference while point_set_2 is data
template <class PointRange1, class PointRange2,
          class NamedParameters1, class NamedParameters2>
#ifdef DOXYGEN_RUNNING
geom_traits::Aff_transformation_3
#else
typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3
#endif
compute_registration_transformation (const PointRange1& point_set_1, const PointRange2& point_set_2,
                                     const NamedParameters1& np1, const NamedParameters2& np2)
{
  using boost::choose_param;
  using boost::get_param;

  namespace PSP = CGAL::Point_set_processing_3;

  // property map types
  typedef typename PSP::GetPointMap<PointRange1, NamedParameters1>::type PointMap1;
  typedef typename PSP::GetPointMap<PointRange2, NamedParameters2>::type PointMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<PointMap1>::value_type,
                                             typename boost::property_traits<PointMap2>::value_type> ::value),
                            "The point type of input ranges must be the same");

  typedef typename PSP::GetNormalMap<PointRange1, NamedParameters1>::type NormalMap1;
  typedef typename PSP::GetNormalMap<PointRange2, NamedParameters2>::type NormalMap2;
  CGAL_static_assertion_msg((boost::is_same< typename boost::property_traits<NormalMap1>::value_type,
                                             typename boost::property_traits<NormalMap2>::value_type> ::value),
                            "The vector type of input ranges must be the same");

  typedef typename PSP::GetK<PointRange1, NamedParameters1>::Kernel Kernel;
  typedef typename Kernel::FT Scalar;
  typedef typename Kernel::Aff_transformation_3 Transformation;

  PointMap1 point_map1 = choose_param(get_param(np1, internal_np::point_map), PointMap1());
  NormalMap1 normal_map1 = choose_param(get_param(np1, internal_np::normal_map), NormalMap1());
  PointMap2 point_map2 = choose_param(get_param(np2, internal_np::point_map), PointMap2());
  NormalMap2 normal_map2 = choose_param(get_param(np2, internal_np::normal_map), NormalMap2());

  // initial transformation
  Transformation initial_transformation
    = choose_param(get_param(np2, internal_np::transformation), Transformation(Identity_transformation()));

  return internal::compute_registration_transformation<Kernel>(point_set_1, point_set_2,
                                                               point_map1, point_map2,
                                                               normal_map1, normal_map2,
                                                               initial_transformation,
                                                               internal::construct_icp<Scalar>(np1, np2));
}

// convenience overloads
template <class PointRange1, class PointRange2,
          class NamedParameters1>
typename CGAL::Point_set_processing_3::GetK<PointRange1, NamedParameters1>
  ::Kernel::Aff_transformation_3
compute_registration_transformation(const PointRange1& point_set_1, const PointRange2& point_set_2,
      const NamedParameters1& np1)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2, np1, params::all_default(point_set_1));
}

template <class PointRange1, class PointRange2>
typename CGAL::Point_set_processing_3::GetK<PointRange1,
          cgal_bgl_named_params<bool, internal_np::all_default_t> >
  ::Kernel::Aff_transformation_3
compute_registration_transformation(const PointRange1& point_set_1, const PointRange2& point_set_2)
{
  namespace params = CGAL::Point_set_processing_3::parameters;
  return compute_registration_transformation(point_set_1, point_set_2,
                                             params::all_default(point_set_1),
                                             params::all_default(point_set_2));
}


} } // end of namespace CGAL::pointmatcher

#endif // CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
