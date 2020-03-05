#include "../include/partition_objects.h"
#include "../include/vars.h"



namespace Skippy {
	Partition_Polyhedron::Partition_Polyhedron(const std::set<Partition_Side> & P)
		: id(++Counters::id_partition_polyhedron)
	{
		for (std::set<Partition_Side>::iterator it_f = P.begin(); it_f != P.end(); it_f++) {
			facets.push_back(*it_f);
			it_f->first->push(this);
		}
	}


	Partition_Polyhedron::~Partition_Polyhedron()
	{
		for (std::list<Partition_Side>::iterator it_f = facets.begin(); it_f != facets.end(); it_f++) {
			it_f->first->pop(this);
		}
	}


	std::list<Partition_Side>::const_iterator Partition_Polyhedron::facets_begin() const
	{
		return facets.cbegin();
	}


	std::list<Partition_Side>::const_iterator Partition_Polyhedron::facets_end() const
	{
		return facets.cend();
	}


	std::list<Partition_Side>::iterator Partition_Polyhedron::facets_begin()
	{
		return facets.begin();
	}


	std::list<Partition_Side>::iterator Partition_Polyhedron::facets_end()
	{
		return facets.end();
	}


	size_t Partition_Polyhedron::facets_size() const
	{
		return facets.size();
	}


	void Partition_Polyhedron::reverse_sides(const int id_p)
	{
		for (std::list<Partition_Side>::iterator it_s = facets.begin() ; it_s != facets.end() ; ++it_s) {
			Partition_Side & s = (*it_s);
			if (s.first->p == id_p) {
				bool eps = s.second;
				s.second = !eps;
			}
		}
	}
}