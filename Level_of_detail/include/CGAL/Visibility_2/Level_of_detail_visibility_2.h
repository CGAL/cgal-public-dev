#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H

// STL includes.
#include <cassert>
#include <map>
#include <vector>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL{

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_visibility_2 {

		public:
			typedef CGAL::Constrained_Delaunay_triangulation_2<KernelTraits> CDT;
			enum class Visibility_label { IN, OUT, UNKNOWN };
			typedef std::map<int, Visibility_label> Visibility_result;

			virtual int compute(const CDT &, const InputContainer &, Visibility_result &) const = 0;
		};

		// This class works only with the xy aligned ground plane that is Plane(0, 0, 1, 0).
		template<class KernelTraits, class InputContainer>
		class Level_of_detail_visibility_from_classification_2 : public Level_of_detail_visibility_2<KernelTraits, InputContainer> {

		public:
			
			typedef Level_of_detail_visibility_2<KernelTraits, InputContainer> Base;
			
			typedef typename Base::Visibility_label  Visibility_label;
			typedef typename Base::Visibility_result Visibility_result;

			typedef KernelTraits   Kernel;
			typedef InputContainer Container;

			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Point_3 Point_3;

			typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel> CDT;
			
			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;

			using Label     = int; 
			using Label_map = typename Container:: template Property_map<Label>; 

			using Visibility_tmp = std::map<int, std::vector<Visibility_label> >;
			using Log = CGAL::LOD::Mylog;

			int compute(const CDT &cdt, const Container &input, Visibility_result &visibility) const {
				
				Log log;
				visibility.clear();

				Visibility_tmp tmp;
				auto number_of_traversed_faces = -1;

				Label_map labels;
				boost::tie(labels, boost::tuples::ignore) = input.template property_map<Label>("label");

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				CGAL::Unique_hash_map<Face_handle, int> F;

				int count = 0;
				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					F[fit] = count++;
					visibility[F[fit]] = Visibility_label::UNKNOWN;
				}

				log.out << "Found faces: " << std::endl;
				for (typename Container::const_iterator it = input.begin(); it != input.end(); ++it) {

					const Point_3 &p = input.point(*it);
					const Point_2 q  = Point_2(p.x(), p.y()); // project p onto xy plane
					
					typename CDT::Locate_type locate_type;
					int locate_index = -1;

					const Face_handle fh = cdt.locate(q, locate_type, locate_index);
					const int face_index = F[fh];

					if (locate_type == CDT::VERTEX ||
						locate_type == CDT::EDGE) {
					
						set_unknown(face_index, tmp);
						continue;
					}

					assert(locate_type != CDT::OUTSIDE_CONVEX_HULL);
					assert(locate_type != CDT::OUTSIDE_AFFINE_HULL);

					assert(locate_type == CDT::FACE);

					const Label label = labels[*it];
					switch (label) {

						case ground:
						set_outside(face_index, tmp);
						break;

						case facade:
						set_inside(face_index, tmp);
						break;

						case roof:
						set_inside(face_index, tmp);
						break;

						case vegetation:
						set_outside(face_index, tmp);
						break;

						default:
						set_unknown(face_index, tmp);
						break;
					}
				}

				postprocess_tmp(tmp, visibility);

				// Remove later.
				for (typename Visibility_result::const_iterator it = visibility.begin(); it != visibility.end(); ++it) {

					const int tmpLabel = static_cast<int>((*it).second);
					std::string labelName = "default";

					if (tmpLabel == 0) labelName = "IN";
					if (tmpLabel == 1) labelName = "OUT";
					if (tmpLabel == 2) labelName = "UNKNOWN";

					log.out << "face index: " << (*it).first << " with label: " << labelName << std::endl;
				}

				number_of_traversed_faces = static_cast<int>(visibility.size());

				log.out << "number of traversed faces: " << number_of_traversed_faces << std::endl;
				log.save("tmp/visibility");

				return number_of_traversed_faces;
			}

			void set_inside(const int face_index, Visibility_tmp &tmp) const {
				tmp[face_index].push_back(Visibility_label::IN);
			}

			void set_outside(const int face_index, Visibility_tmp &tmp) const {
				tmp[face_index].push_back(Visibility_label::OUT);
			}

			void set_unknown(const int face_index, Visibility_tmp &tmp) const {
				tmp[face_index].push_back(Visibility_label::UNKNOWN);
			}

			void postprocess_tmp(const Visibility_tmp &tmp, Visibility_result &visibility) const {

				for (typename Visibility_tmp::const_iterator it = tmp.begin(); it != tmp.end(); ++it) {

					int in_count      = 0;
					int out_count     = 0;
					int unknown_count = 0;

					for (size_t i = 0; i < (*it).second.size(); ++i) {
						switch ((*it).second[i]) {

							case Visibility_label::IN:
								++in_count;
								break;

							case Visibility_label::OUT:
								++out_count;
								break;

							case Visibility_label::UNKNOWN:
								++unknown_count;
								break;

							default:
								break;
						}
					}

					assert(in_count + out_count + unknown_count != 0);

					if (in_count >= out_count && in_count >= unknown_count) {
						visibility[(*it).first] = Visibility_label::IN;
						continue;
					}

					if (out_count >= in_count && out_count >= unknown_count) {
						visibility[(*it).first] = Visibility_label::OUT;
						continue;
					}
					
					if (unknown_count >= in_count && unknown_count >= out_count) {
						visibility[(*it).first] = Visibility_label::UNKNOWN;
						continue;
					}

					const bool continue_failure = false;
					assert(continue_failure);
				}
			}
		};		
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H