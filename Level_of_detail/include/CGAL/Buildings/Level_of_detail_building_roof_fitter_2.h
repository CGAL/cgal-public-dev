#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Main class.
		template<class KernelTraits, class CDTInput, class ContainerInput, class FittingStrategyInput>
		class Level_of_detail_building_roof_fitter_2 {

		public:
			typedef KernelTraits    	 Kernel;
			typedef CDTInput        	 CDT;
			typedef ContainerInput  	 Container;
			typedef FittingStrategyInput FittingStrategy;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_3 Point_3;

			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;

			// Extra.
			using Point_index = typename Container::Index;
			using Face_points_map = std::map<Face_handle, std::vector<Point_index> >;

			using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
			using Buildings = std::map<int, Building>;

			using Building_iterator = typename Buildings::iterator;

			using Label     = int; 
			using Label_map = typename Container:: template Property_map<Label>;

			Level_of_detail_building_roof_fitter_2() { }

			void fit_roof_heights(const CDT &, const Container &input, const Face_points_map &fp_map, Buildings &buildings) {
				
				set_input(input);

				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
					
					m_fitter.clear();
					Building &building = (*bit).second;

					const FT height = fit_height_for_building(input, fp_map, building);
					add_height_to_building(height, building);
				}
			}

		private:
			Label_map m_labels;
			FittingStrategy m_fitter;

			void set_input(const Container &input) {
				boost::tie(m_labels, boost::tuples::ignore) = input.template property_map<Label>("label");
			}

			FT fit_height_for_building(const Container &input, const Face_points_map &fp_map, const Building &building) {

				const size_t num_faces = building.faces.size();
				for (size_t i = 0; i < num_faces; ++i) add_face_heights(input, fp_map, building.faces[i]);

				return m_fitter.get_result();
			}

			void add_face_heights(const Container &input, const Face_points_map &fp_map, const Face_handle &fh) {

				// assert(is_valid_face(fp_map, fh));   // may be slow, if I add all fhs to the fp_map, I do not need to do it at all
				if (!is_valid_face(fp_map, fh)) return; // if fh is not in the fp_map, it is probably because this fh does not have any associated points

				const size_t num_points = fp_map.at(fh).size();
				for (size_t i = 0; i < num_points; ++i) {
					
					const Point_index point_index = fp_map.at(fh)[i];
					const Label point_label = m_labels[point_index];

					if (is_valid_point(point_label)) {

						const Point_3 &p = input.point(point_index);						
						const FT height = p.z();

						m_fitter.add_height(height);
					}
				}
			}

			bool is_valid_face(const Face_points_map &fp_map, const Face_handle &fh) {
				return fp_map.count(fh);
			}

			bool is_valid_point(const Label point_label) {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				switch (point_label) {

					case ground: 	 return false;
					case facade: 	 return false;
					case roof: 		 return true;
					case vegetation: return false;

					default: assert(!"Wrong label!");
				}
			}

			void add_height_to_building(const FT height, Building &building) {
				building.height = height;
			}
		};


		// Fitter base.
		template<class KernelTraits>
		class Level_of_detail_height_fitter_base {
		
		public:
			typedef KernelTraits        Kernel;
			typedef typename Kernel::FT FT;

			virtual void add_height(const FT value) = 0;
			virtual FT   get_result() = 0;
			virtual void clear() = 0;

			virtual ~Level_of_detail_height_fitter_base() { }
		};


		// Min height fitter.
		template<class KernelTraits>
		class Level_of_detail_min_height_fitter : public Level_of_detail_height_fitter_base<KernelTraits> {

		public:
			typedef Level_of_detail_height_fitter_base<KernelTraits> Base;

			typedef typename Base::Kernel Kernel;
			typedef typename Base::FT     FT;

			Level_of_detail_min_height_fitter() : m_big_value(FT(1000000)), m_min_height(m_big_value) { }

			void add_height(const FT value) override {

				assert(value > FT(0));
				m_min_height = CGAL::min(m_min_height, value);
			}

			FT get_result() override {
				return m_min_height;
			}

			void clear() override {
				m_min_height = m_big_value;
			}

		private:
			FT m_big_value;
			FT m_min_height;
		};


		// Average height fitter.
		template<class KernelTraits>
		class Level_of_detail_avg_height_fitter : public Level_of_detail_height_fitter_base<KernelTraits> {

		public:
			typedef Level_of_detail_height_fitter_base<KernelTraits> Base;

			typedef typename Base::Kernel Kernel;
			typedef typename Base::FT     FT;

			Level_of_detail_avg_height_fitter() : m_num_values(FT(0)), m_sum_height(FT(0)) { }

			void add_height(const FT value) override {

				assert(value > FT(0));
				
				m_sum_height += value;
				m_num_values += FT(1);
			}

			FT get_result() override {
				return m_sum_height / m_num_values;
			}

			void clear() override {

				m_sum_height = FT(0);
				m_num_values = FT(0);
			}

		private:
			FT m_num_values;
			FT m_sum_height;
		};


		// Max height fitter.
		template<class KernelTraits>
		class Level_of_detail_max_height_fitter : public Level_of_detail_height_fitter_base<KernelTraits> {

		public:
			typedef Level_of_detail_height_fitter_base<KernelTraits> Base;

			typedef typename Base::Kernel Kernel;
			typedef typename Base::FT     FT;

			Level_of_detail_max_height_fitter() : m_negative_value(-FT(1)), m_max_height(m_negative_value) { }

			void add_height(const FT value) override {

				assert(value > FT(0));
				m_max_height = CGAL::max(m_max_height, value);
			}

			FT get_result() override {
				return m_max_height;
			}

			void clear() override {
				m_max_height = m_negative_value;
			}

		private:
			FT m_negative_value;
			FT m_max_height;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FITTER_2_H