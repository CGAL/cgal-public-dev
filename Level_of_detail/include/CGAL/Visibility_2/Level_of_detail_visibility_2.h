#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H

// STL includes.
#include <utility>
#include <cassert>
#include <vector>
#include <map>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/Barycentric_coordinates_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_2 {

		public:
			typedef KernelTraits 		Kernel;
			typedef typename Kernel::FT FT;

			typedef ContainerInput Container;
			typedef CDTInput 	   CDT;
			
			virtual int compute(const Container &, CDT &) const = 0;
		};

		// This class works only with the xy aligned ground plane that is Plane(0, 0, 1, 0).
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_from_classification_2 : public Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			typedef Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> Base;

			typedef typename Base::Kernel 	 Kernel;
			typedef typename Base::FT     	 FT;
			typedef typename Kernel::Point_2 Point_2;

			typedef typename Base::Container  Container;			
			typedef typename Base::CDT 		  CDT;
			
			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;

			typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

			typedef int Label;
			typedef typename std::pair<Point_2, Label> 							Point_with_label;
			typedef typename CGAL::First_of_pair_property_map<Point_with_label> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					  Search_traits_2;
			typedef CGAL::Search_traits_adapter<Point_with_label, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 						  Neighbor_search;
			typedef typename Neighbor_search::Tree 											  Tree;

			using Visibility = std::map<Face_handle, std::vector<Visibility_label> >;
			using Log = CGAL::LOD::Mylog;

			using Point_iterator = typename Container::const_iterator;

			Level_of_detail_visibility_from_classification_2() : 
			m_approach(Visibility_approach::POINT_BASED), 
			m_method(Visibility_method::CLASSIFICATION), 
			m_sampler(Visibility_sampler::REGULAR),
			m_num_samples(1),
			m_k(6) { }

			int compute(const Container &input, CDT &cdt) const {

				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) fit->info().in = 0.5;

				switch(m_approach) {

					case Visibility_approach::POINT_BASED:
					compute_point_based_visibility(input, cdt);
					break;

					case Visibility_approach::FACE_BASED:
					compute_face_based_approach(input, cdt);
					break;

					default:
					assert("Wrong approach!");
					break;
				}

				// Remove later.
				Log log;
				log.out << "Visibility labels: " << std::endl;

				int count = 0;
				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++count) {

					const FT result = fit->info().in;
					std::string labelName = "default";

					const FT half = FT(1) / FT(2);

					if (result > half)  labelName = "IN";
					if (result < half)  labelName = "OUT";
					if (result == half) labelName = "UNKNOWN";

					log.out << "face index: " << count << " with label: " << labelName << " and visibility: " << result << std::endl;
				}
				log.save("tmp/visibility");

				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			const Visibility_approach m_approach;
			const Visibility_method   m_method;
			const Visibility_sampler  m_sampler;

			const size_t m_num_samples;
			const size_t m_k; 			// change it to autodetection later!

			void compute_point_based_visibility(const Container &input, CDT &cdt) const {

				Visibility visibility;
				for (typename Container::const_iterator it = input.begin(); it != input.end(); ++it) {

					const Point_2 &p = (*it).first;
					typename CDT::Locate_type locate_type;
					
					int locate_index = -1;
					const Face_handle face_handle = cdt.locate(p, locate_type, locate_index);

					if (locate_type == CDT::VERTEX ||
						locate_type == CDT::EDGE /*||
						on_the_border(cdt, fh, p) */) {
					
						// Improve this part of the code if possible!
						set_unknown(face_handle, visibility);
						continue;
					}

					assert(locate_type != CDT::OUTSIDE_CONVEX_HULL);
					assert(locate_type != CDT::OUTSIDE_AFFINE_HULL);

					assert(locate_type == CDT::FACE);

					switch(m_method) {

						case Visibility_method::CLASSIFICATION: {

							const Label point_label = (*it).second;
							estimate_with_classification(face_handle, point_label, visibility);
							break;
						}

						default:
						assert("Wrong visibility method!");
						break;
					}
				}
				postprocess(visibility);
			}

			void estimate_with_classification(const Face_handle face_handle, const Label point_label, Visibility &visibility) const {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				switch (point_label) {

					case ground:
					set_outside(face_handle, visibility);
					break;

					case facade:
					set_unknown(face_handle, visibility);
					break;

					case roof:
					set_inside(face_handle, visibility);
					break;

					case vegetation:
					set_outside(face_handle, visibility);
					break;

					default:
					assert("Classification label is missing!");
					break;
				}
			}

			void compute_face_based_approach(const Container &input, CDT &) const {

				// Remove it later and use property maps instead.
				/*
				std::vector<Point_2> points(input.number_of_points());
				for (Point_iterator it = input.begin(); it != input.end(); ++it) {
					
					const int index = static_cast<int>(*it);
					const Point_3 &p = input.point(index);

					points[index] = Point_2(p.x(), p.y()); // project onto xy plane
				} */

				// Create a tree.
				Tree tree(input.begin(), input.end());

				Point_2 query(FT(0), FT(0));
				Neighbor_search search(tree, query, m_k);

				// Iterate over all faces.
				/*
				std::vector<Point_2> samples(m_num_samples);
				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

					generate_samples(cdt, fit, samples);
					assert(samples.size() == m_num_samples);

					std::vector<FT> visibility(m_num_samples);
					compute_visibilities(tree, points);
				} */
			}

			void generate_samples(std::vector<Point_2> &samples) const {
				
				switch(m_sampler) {

					case Visibility_sampler::REGULAR:
						generate_regular_samples(samples);
						break;

					case Visibility_sampler::RANDOM:
						break;

					default:
						assert("Wrong sampler!");
						break;
				}
			}

			void generate_regular_samples(const CDT &, const Face_handle &, std::vector<Point_2> &samples) const {

				samples[0] = Point_2(0, 0);
			}

			void set_inside(const Face_handle face_handle, Visibility &visibility) const {
				visibility[face_handle].push_back(Visibility_label::IN);
			}

			void set_outside(const Face_handle face_handle, Visibility &visibility) const {
				visibility[face_handle].push_back(Visibility_label::OUT);
			}

			void set_unknown(const Face_handle face_handle, Visibility &visibility) const {
				visibility[face_handle].push_back(Visibility_label::UNKNOWN);
			}

			void postprocess(Visibility &visibility) const {

				for (typename Visibility::const_iterator it = visibility.begin(); it != visibility.end(); ++it) {
					
					FT inside  = FT(0);
					FT outside = FT(0);

					const FT half = FT(1) / FT(2);
					for (size_t i = 0; i < (*it).second.size(); ++i) {

						switch ((*it).second[i]) {

							case Visibility_label::IN:
								inside  += FT(1);
								break;

							case Visibility_label::OUT:
								outside += FT(1);
								break;

							case Visibility_label::UNKNOWN: {

								inside  += half;
								outside += half;

								break;
							}

							default:
								assert("Should never get here!");
								break;
						}
					}

					const FT sum = FT(inside + outside); 
					if (sum == FT(0)) {

						(*it).first->info().in = half;
						continue;	
					}
					(*it).first->info().in = inside / sum;
				}
			}

			bool on_the_border(const CDT &cdt, const Face_handle fh, const Point_2 &p) const {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				Triangle_coordinates tric(a, b, c);
				std::vector<FT> bc;
				tric(p, bc);

				const FT tol = 1.0e-6;

				if (bc[0] < tol || bc[1] < tol || bc[2] < tol) return true;
				return false;
			}
		};		
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H