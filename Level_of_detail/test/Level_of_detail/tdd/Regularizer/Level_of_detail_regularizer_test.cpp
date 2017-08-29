// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Regularizer/Level_of_detail_regularizer.h>

using namespace testing;

class LOD_RegularizerTest: public Test {
	
public:
	using FT = double;
	
	using Traits = CGAL::Simple_cartesian<FT>;
	
	using Point  = Traits::Point_3;
	using Normal = Traits::Vector_3;
	using Plane  = Traits::Plane_3;

	using Planes         = std::map<int, std::vector<int> >;
	using Container      = CGAL::Point_set_3<Point>;
	using LodRegularizer = CGAL::LOD::Level_of_detail_vertical_regularizer<Traits, Container, Planes>;

	using Plane_map = Container:: template Property_map<Plane>;
	using Iterator  = Container::iterator;

	LodRegularizer lodRegularizer;

	void set_input_properties(Container &input, Plane_map &planes) const {
		
		Point zero(0, 0, 0);
		Plane default_plane(zero, zero, zero);

		input.add_normal_map();

		auto success = false;
		boost::tie(planes, success)  = input. template add_property_map<Plane>("plane", default_plane);
		assert(success);
	}

	void create_test_input(Container &input, Planes &indices, Plane_map &planes, Plane &ground) const {

		set_input_properties(input, planes);

		Iterator 
		it = input.insert(Point(0.2, 0.3, 0.22), Normal(0.0, 0.22, 0.0)); planes[*it] = Plane(0.0, 0.22, 0.0, 0.0);
		it = input.insert(Point(0.4, 0.6, 0.22), Normal(0.0, 0.22, 0.0)); planes[*it] = Plane(0.0, 0.22, 0.0, 0.0);
		it = input.insert(Point(0.8, 0.1, 0.22), Normal(0.0, 0.22, 0.0)); planes[*it] = Plane(0.0, 0.22, 0.0, 0.0);

		// Plane to be regularized!
		it = input.insert(Point(1.00, 0.2, 0.9), Normal(0.36, 0.0, 0.03)); planes[*it] = Plane(0.36, 0.0, 0.03, 0.387);
		it = input.insert(Point(1.05, 0.5, 0.3), Normal(0.36, 0.0, 0.03)); planes[*it] = Plane(0.36, 0.0, 0.03, 0.387);
		it = input.insert(Point(1.00, 0.8, 0.9), Normal(0.36, 0.0, 0.03)); planes[*it] = Plane(0.36, 0.0, 0.03, 0.387);

		indices[0].push_back(0); indices[0].push_back(1); indices[0].push_back(2);
		indices[1].push_back(3); indices[1].push_back(4); indices[1].push_back(5);

		ground = Plane(0.0, 0.0, 1.0, 0.0);
	}
};

TEST_F(LOD_RegularizerTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_RegularizerTest, UnregularizedPlaneIsPreserved) {

	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	lodRegularizer.regularize(indices, input, ground);
	
	ASSERT_THAT(input.point(0) , Eq(Point(0.2, 0.3, 0.22)));
	ASSERT_THAT(input.normal(1), Eq(Normal(0.0, 0.22, 0.0)));
	ASSERT_THAT(planes[2]      , Eq(Plane(0.0, 0.22, 0.0, 0.0)));
}

TEST_F(LOD_RegularizerTest, ReturnsOneRegularizedPlane) {

	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	const auto number_of_regularized_planes = lodRegularizer.regularize(indices, input, ground);

	ASSERT_THAT(static_cast<int>(number_of_regularized_planes), Eq(1));
}

TEST_F(LOD_RegularizerTest, RegularizedPointsHaveEqualNormals) {

	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	lodRegularizer.regularize(indices, input, ground);

	ASSERT_THAT(input.normal(3), Eq(input.normal(4)));
	ASSERT_THAT(input.normal(4), Eq(input.normal(5)));
}

TEST_F(LOD_RegularizerTest, RegularizedPointsHaveEqualPlanes) {

	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	lodRegularizer.regularize(indices, input, ground);

	boost::tie(planes, boost::tuples::ignore) = input. template property_map<Plane>("plane");

	ASSERT_THAT(planes[3], Eq(planes[4]));
	ASSERT_THAT(planes[4], Eq(planes[5]));
}

TEST_F(LOD_RegularizerTest, ChangesPointsNormalsAndPlanesAfterRegularization) {

	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	lodRegularizer.regularize(indices, input, ground);
	
	const auto point_diff  = CGAL::squared_distance(input.point(3) , Point(1.07129, 0.2, 0.813846));
	const auto normal_diff = (input.normal(4) - Normal(0.361248, 0.0, 0.0)).squared_length();
	
	const auto a = CGAL::abs(planes[5].a() - 0.361248);
	const auto b = CGAL::abs(planes[5].b() - 0.0);
	const auto c = CGAL::abs(planes[5].c() - 0.0);
	const auto d = CGAL::abs(planes[5].d() + 0.387);

	const auto aa = a * a;
	const auto bb = b * b;
	const auto cc = c * c;
	const auto dd = d * d;

	const auto plane_diff = aa + bb + cc + dd;

	const auto eps = 1.0e-6;

	ASSERT_LT(point_diff , eps);
	ASSERT_LT(normal_diff, eps);
	ASSERT_LT(plane_diff , eps);
}

TEST_F(LOD_RegularizerTest, UsesAverageNormalAsPlaneNormal) {
	
	Plane ground;
	Plane_map planes; Planes indices; Container input;

	create_test_input(input, indices, planes, ground);
	lodRegularizer.regularize(indices, input, ground);

	const auto average_normal = Normal(0.361248, 0.0, 0.0);
	const auto normal_diff = (input.normal(5) - average_normal).squared_length();

	const auto eps = 1.0e-6;

	ASSERT_LT(normal_diff, eps);
}
