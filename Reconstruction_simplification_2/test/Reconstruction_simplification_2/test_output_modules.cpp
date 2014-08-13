// test_output_modules.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include <CGAL/List_output.h>
#include <CGAL/Index_output.h>
#include <CGAL/Tds_output.h>

#include<fstream>
#include<iostream>
#include <string>
#include <cassert>
#include <iterator>
#include <utility>      // std::pair


#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>

#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                     	                    Point;
typedef K::Segment_2                 						Segment;

typedef K::FT                                         		FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair>  PointPMap;
typedef CGAL::Second_of_pair_property_map <PointMassPair> MassPMap;


typedef CGAL::Reconstruction_simplification_2<K, PointPMap, MassPMap> Rs_2;

typedef Rs_2::Vertex Vertex;

typedef Rs_2::Reconstruction_edge_2 R_edge_2;

typedef CGAL::Reconstruction_triangulation_2<K> Rt_2;

typedef Rt_2::Finite_edges_iterator Finite_edges_iterator;
typedef Rt_2::Vertex_iterator Vertex_iterator;

typedef Rt_2::Edge Edge;


void test_list_output(Rs_2& rs2);
void test_index_output(Rs_2& rs2);
void test_tds_output(Rs_2& rs2);


void print_edge(Edge edge) {
	int i = edge.second;
	Point a = edge.first->vertex((i+1)%3)->point();
	Point b = edge.first->vertex((i+2)%3)->point();
	std::cout << a << " , " << b <<  std::endl;
}

int main ()
{

	PointMassList points;
	//use the stair example for testing
	load_xy_file<PointMassList, Point>("data/stair-noise00.xy", points);

	PointPMap point_pmap;
	MassPMap  mass_pmap;

	Rs_2 rs2(points.begin(), points.end(), point_pmap, mass_pmap);

	rs2.reconstruct(100); //100 steps

	test_list_output(rs2);
	test_tds_output(rs2);
	test_index_output(rs2);
}

void test_list_output(Rs_2& rs2) {

	std::cout <<"(-------------List OUTPUT---------- )" << std::endl;

	std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	typedef std::back_insert_iterator<std::vector<Point> > Point_it;
	typedef std::back_insert_iterator<std::vector<Segment> >  Edge_it;

	Point_it point_it(isolated_points);
	Edge_it  edge_it(edges);

	CGAL::List_output<K, Point_it, Edge_it> list_output(point_it, edge_it);

    rs2.extract_solid_elements(list_output);

    int vertex_count = 0;
	for (std::vector<Point>::iterator it = isolated_points.begin();
			it != isolated_points.end(); it++) {
		vertex_count++;
  		std::cout  <<  *it << std::endl;
	}
	assert(vertex_count == 18);

	int edge_count = 0;
	for (std::vector<Segment>::iterator it = edges.begin();
			it != edges.end(); it++) {
  		std::cout << *it << std::endl;
  		edge_count++;
    }
	assert(edge_count == 31);


}

void test_tds_output(Rs_2& rs2) {

	std::cout <<"(-------------Tds OUTPUT---------- )" << std::endl;

	CGAL::Tds_output<K> tds_output;
	rs2.extract_solid_elements(tds_output);
	Rt_2 rt2;
	tds_output.extract_reconstruction_tds(rt2);

	int vertex_count = 0;
	for (Vertex_iterator vi = rt2.vertices_begin();
					  vi != rt2.vertices_end(); ++vi) {

		FT relevance = (*vi).get_relevance();
		if (relevance > 0) {
			std::cout  <<  *vi << std::endl;
			vertex_count++;
		}
	}
	std::cout <<"vertex_count " << vertex_count << std::endl;
	assert(vertex_count == 18);

	int edge_count = 0;
	for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei) {
		FT relevance = (*ei).first->relevance((*ei).second);
		if (relevance > 0) {
			print_edge(*ei);

			std::cout << "mass "   << tds_output.get_mass(*ei) << std::endl;
			std::cout << "cost "   << tds_output.get_cost(*ei) << std::endl;
			std::cout << "length " << tds_output.get_length(*ei) << std::endl;

			edge_count++;
		}
	}
	std::cout <<"edge_count " << edge_count << std::endl;
	assert(edge_count == 31);

}
void test_index_output(Rs_2& rs2) {

	std::cout <<"(-------------Off OUTPUT---------- )" << std::endl;

	CGAL::Index_output<K> index_output;

    rs2.extract_solid_elements(index_output);

    //print
    index_output.get_os_output(std::cout);

    //test cardinalities
    std::ostringstream buffer;
    index_output.get_os_output(buffer);

    std::stringstream stream(buffer.str());
	std::vector<std::string> res;
	while (1){
		std::string line;
		std::getline(stream,line);
		if (!stream.good())
			break;
		res.push_back(line);
	}

	assert(res.size() == 110);

	assert(res.front() == "OFF 60 0 31");

	for (int i = 61; i < 79; i++) {
		assert(res[i].substr(0,2) == "1 ");
	}

	for (int i = 79; i < 110; i++) {
		assert(res[i].substr(0,2) == "2 ");
	}
}
