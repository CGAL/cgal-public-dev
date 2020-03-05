#include "../include/polygon_vertex_octree.h"
#include "../include/octree_base.h"
#include "../include/octree_base_vertex.h"



namespace Skippy {

	int Polygon_Vertex_Octree::id_vertex = -1;

	Polygon_Vertex_Octree::Polygon_Vertex_Octree(const std::vector<double> & dims)
		: Octree_Base(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], nullptr)
	{

	}


	Polygon_Vertex_Octree::~Polygon_Vertex_Octree()
	{
		
	}



	int Polygon_Vertex_Octree::get_identifier(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M)
	{
		// Search for a point defined by M and approximated by M.
		// If we can find such a point, returns its identifer.
		// Otherwise, adds an element to the octree.

		std::list<Octree_Base_Vertex*> R;
		search(M, hint_M, R);

		for (auto it_v = R.begin() ; it_v != R.end() ; ++it_v) {
			if ((*it_v)->M == M) {
				Polygon_Vertex_Octree_Data* v = dynamic_cast<Polygon_Vertex_Octree_Data*>(*it_v);
				return v->id;
			}
		}

		Polygon_Vertex_Octree_Data* v = new Polygon_Vertex_Octree_Data(M, hint_M, ++id_vertex);
		add(v);

		return v->id;
	}


	void Polygon_Vertex_Octree::get_all_sorted_vertices(std::vector<Polygon_Vertex_Octree_Data*> & R) const
	{
		std::list<Octree_Base_Vertex*> R_0;
		get_all_vertices(R_0);

		R.reserve(R_0.size());
		for (std::list<Octree_Base_Vertex*>::const_iterator it_v = R_0.begin() ; it_v != R_0.end() ; ++it_v) {
			R.push_back(dynamic_cast<Polygon_Vertex_Octree_Data*>(*it_v));
		}

		struct _Vertex_Comparator {
			bool operator() (const Polygon_Vertex_Octree_Data* v1, const Polygon_Vertex_Octree_Data* v2) {
				return v1->id < v2->id;
			}
		} Vertex_Comparator;

		std::sort(R.begin(), R.end(), Vertex_Comparator);
	}
}