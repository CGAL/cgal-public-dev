
namespace CGAL {

	template <typename Kernel>
	Hypothesis<Kernel>::~Hypothesis() {
		for (std::size_t i = 0; i < supporting_planes_.size(); ++i)
			delete supporting_planes_[i];
		supporting_planes_.clear();
	}


	template <typename Kernel>
	void Hypothesis<Kernel>::generate(Mesh& mesh) {
		refine_planes();

		Mesh bbox_mesh;
		construct_bbox_mesh(bbox_mesh);

		construct_proxy_mesh(bbox_mesh, mesh);

		pairwise_intersection(mesh);

		typedef std::set<Plane*>	Plane_set;
		// store where (i.e., the two planes) an edge is computed from. 
		Mesh::Property_map< Edge_descriptor, Plane_set >		edge_source_planes_;

		// store where (i.e., the plane triplet) a vertex is computed from.
		typename Mesh::Property_map< Vertex_descriptor, Plane_set >		vertex_source_planes_;

	}

	namespace details {
		template <typename Planar_segment>
		class SegmentSizeCmpInc
		{
		public:
			SegmentSizeCmpInc() {}
			bool operator()(Planar_segment* s0, Planar_segment* s1) const {
				return s0->size() < s1->size();
			}
		};

		template <typename FT, typename Vector>
		void normalize(Vector& v) {
			FT s = std::sqrt(v.squared_length());
			if (s > 1e-30)
				s = FT(1) / s;
			v *= s;
		}

		template <typename Point_set, typename Planar_segment, typename Plane, typename FT>
		std::size_t num_points_on_plane(Point_set* point_set, const Planar_segment* s, const Plane& plane, FT dist_threshold) {
			typedef typename Point_set::FT			FT;
			typedef typename Point_set::Point		Point;
			typedef typename Point_set::Vector		Vector;

			std::size_t count = 0;
			const Point_set::Point_map& points = point_set->point_map();
			for (std::size_t i = 0; i < s->size(); ++i) {
				std::size_t idx = s->at(i);
				const Point& p = points[idx];

				FT sdist = CGAL::squared_distance(plane, p);
				FT dist = std::sqrt(sdist);
				if (dist < dist_threshold)
					++ count;
			}
			return count;
		}

		template <typename Point_set, typename Planar_segment, typename Plane, typename FT>
		void merge(Point_set* point_set, Planar_segment* g1, Planar_segment* g2, FT max_dist) {
			std::vector< Planar_segment* >& groups = point_set->planar_segments();

			std::vector<std::size_t> points_indices;
			points_indices.insert(points_indices.end(), g1->begin(), g1->end());
			points_indices.insert(points_indices.end(), g2->begin(), g2->end());

			Planar_segment* g = new Planar_segment;
			g->insert(g->end(), points_indices.begin(), points_indices.end());
			g->set_point_set(point_set);
			g->fit_plane();
			groups.push_back(g);

			std::vector< Planar_segment* >::iterator pos = std::find(groups.begin(), groups.end(), g1);
			if (pos != groups.end()) {
				Planar_segment* tmp = *pos;
				groups.erase(pos);
				delete tmp;
			}

			pos = std::find(groups.begin(), groups.end(), g2);
			if (pos != groups.end()) {
				Planar_segment* tmp = *pos;
				groups.erase(pos);
				delete tmp;
			}
		}

	}

	template <typename Kernel>
	void Hypothesis<Kernel>::refine_planes() {
		Point_set* pset = const_cast<Point_set*>(point_set_);
		std::vector< Planar_segment* >& segments = pset->planar_segments();
		const Point_set::Point_map& points = pset->point_map();

		FT avg_max_dist = 0;
		for (std::size_t i = 0; i < segments.size(); ++i) {
			Planar_segment* s = segments[i];
			s->fit_plane();
			const Plane& plane = s->supporting_plane();

			FT g_max_dist = -FLT_MAX;
			for (std::size_t j = 0; j < s->size(); ++j) {
				std::size_t idx = s->at(j);
				const Point& p = points[idx];
				FT sdist = CGAL::squared_distance(plane, p);
				g_max_dist = std::max(g_max_dist, std::sqrt(sdist));
			}

			avg_max_dist += g_max_dist;
		}
		avg_max_dist /= segments.size();
		avg_max_dist /= 2.0f;

		FT theta = 10.0f;				// in degree
		theta = static_cast<FT>(CGAL_PI * theta / 180.0f);	// in radian
		bool merged = false;
		do
		{
			merged = false;
			std::sort(segments.begin(), segments.end(), details::SegmentSizeCmpInc<Planar_segment>());

			for (std::size_t i = 0; i < segments.size(); ++i) {
				Planar_segment* g1 = segments[i];
				const Plane& plane1 = g1->supporting_plane();
				Vector n1 = plane1.orthogonal_vector();	details::normalize<FT, Vector>(n1);
				FT num_threshold = g1->size() / 5.0f;
				for (std::size_t j = i + 1; j < segments.size(); ++j) {
					Planar_segment* g2 = segments[j];
					const Plane& plane2 = g2->supporting_plane();
					Vector n2 = plane2.orthogonal_vector();	details::normalize<FT, Vector>(n2);
					if (std::abs(n1 * n2) > std::cos(theta)) {
						std::size_t set1on2 = details::num_points_on_plane<Point_set, Planar_segment, Plane, FT>(pset, g1, plane2, avg_max_dist);
						std::size_t set2on1 = details::num_points_on_plane<Point_set, Planar_segment, Plane, FT>(pset, g2, plane1, avg_max_dist);
						if (set1on2 > num_threshold || set2on1 > num_threshold) {
							details::merge<Point_set, Planar_segment, Plane, FT>(pset, g1, g2, avg_max_dist);
							merged = true;
							break;
						}
					}
				}
				if (merged)
					break;
			}
		} while (merged);

		std::sort(segments.begin(), segments.end(), details::SegmentSizeCmpInc<Planar_segment>());
	}


	template <typename Kernel>
	void Hypothesis<Kernel>::construct_bbox_mesh(Mesh& mesh) {

	}


	template <typename Kernel>
	void Hypothesis<Kernel>::construct_proxy_mesh(const Mesh& bbox_mesh, Mesh& mesh) {

	}


	template <typename Kernel>
	void Hypothesis<Kernel>::pairwise_intersection(Mesh& mesh) {

	}

}