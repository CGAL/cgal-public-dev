#pragma once
#include "defs.h"
#include "propagation.h"



namespace Skippy
{
	class Kinetic_Propagation_Simple : public Kinetic_Propagation
	{
	public:
		KINETIC_PARTITION_API Kinetic_Propagation_Simple();

		KINETIC_PARTITION_API Kinetic_Propagation_Simple(const int N);

		KINETIC_PARTITION_API Kinetic_Propagation_Simple(int argc, char *argv[]);

		KINETIC_PARTITION_API Kinetic_Propagation_Simple(const std::string & filename, Preprocess process = NONE);

		KINETIC_PARTITION_API Kinetic_Propagation_Simple(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process = NONE);

		KINETIC_PARTITION_API Kinetic_Propagation_Simple(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process = NONE);

		KINETIC_PARTITION_API virtual ~Kinetic_Propagation_Simple();

		KINETIC_PARTITION_API void run();

		void delete_unique_kinetic_data_structure();

	protected:
		void init_unique_kinetic_data_structure();

		void build_polygons() const;

		void build_partition();
	};
}