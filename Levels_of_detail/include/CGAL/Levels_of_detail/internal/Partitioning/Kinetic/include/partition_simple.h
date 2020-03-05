#pragma once
#include "partition.h"


namespace Skippy 
{
	class Partition_Simple : public Partition 
	{
	public:
		Partition_Simple(const std::vector<CGAL_Plane> & plane_defs);

		~Partition_Simple();

		void build();
	};
}