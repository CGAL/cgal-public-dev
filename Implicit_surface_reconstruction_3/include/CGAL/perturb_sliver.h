#ifndef CGAL_IMPLICIT_RECONSTRUCTION_PERTURB_SLIVER_H
#define CGAL_IMPLICIT_RECONSTRUCTION_PERTURB_SLIVER_H

#include <unordered_map>
#include <vector>

#include <CGAL/Modifiable_priority_queue.h>

#include <boost/shared_ptr.hpp>

namespace CGAL{

    template<class Gt,
             class Tr>
    class Sliver_perturbation_removal
    {
        typedef Gt Geom_traits;
        typedef Tr Triangulation;
        
        // Geometric types
        typedef typename Geom_traits::Point_3 Point; ///< point type.
        typedef typename Geom_traits::FT FT;
        typedef typename Geom_traits::Vector_3 Vector;

        // Triangulation types
        typedef typename Triangulation::Edge                         Edge;
        typedef typename Triangulation::Vertex_handle                Vertex_handle;
        typedef typename Triangulation::Cell_handle                  Cell_handle;
        typedef typename Triangulation::Finite_cells_iterator        Finite_cells_iterator;
    
        typedef typename std::pair<std::vector<Cell_handle>, FT>     Sliver_info;

    private:

        // class member
        boost::shared_ptr<Triangulation> & m_tr;
        bool is_octree;

        enum PerturbationMethod { NULL_METHOD, OCTREE, RADIUS, VOLUME, RANDOM };

        template< typename FT
            , typename Vertex_handle
            , typename Point_3
            , typename PerturbationMethod>
            class PVertex_
        {
        public:
            typedef PVertex_<FT,
                Vertex_handle,
                Point_3,
                PerturbationMethod> Self;

            /// Constructor
            PVertex_()
                : vertex_handle_()
                , incident_sliver_nb_(0)
                , min_value_((std::numeric_limits<double>::max)())
                , try_nb_(0)
                , p_perturbation_(NULL_METHOD)
                , id_()
                , is_interior_()
            { }

            PVertex_(const Vertex_handle& vh, size_t id)
                : vertex_handle_(vh)
                , incident_sliver_nb_(0)
                , min_value_((std::numeric_limits<double>::max)())
                , try_nb_(0)
                , p_perturbation_(NULL_METHOD)
                , id_(id)
                , is_interior_(this->has_finite_voronoi_cell(vh))
            { }

            /// Associated vertex
            const Vertex_handle& vertex() const { return vertex_handle_; }
            void set_vertex(const Vertex_handle& vh) { vertex_handle_ = vh; }

            /// Incident slivers number
            unsigned int sliver_nb() const { return incident_sliver_nb_; }
            void set_sliver_nb(const unsigned int n) { incident_sliver_nb_ = n; }

            /// Current perturbation
            const PerturbationMethod perturbation() const { return p_perturbation_; }
            void set_perturbation(const PerturbationMethod p) { p_perturbation_ = p; }
            void next_perturbation() {
                if (p_perturbation_ == OCTREE)
                    p_perturbation_ = RADIUS;
                else if (p_perturbation_ == VOLUME)
                    p_perturbation_ = RANDOM;
                else if (p_perturbation_ == VOLUME)
                    p_perturbation_ = RANDOM;
            }

            /// Is perturbable
            bool is_perturbable() const
            {  
                return (!(NULL_METHOD == perturbation()) && (sliver_nb() != 0) );
            }

            /// Min sliver value
            const FT& min_value() const { return min_value_; }
            void set_min_value(const FT& min_value){ min_value_ = min_value; }

            /// Try nb
            const unsigned int& try_nb() const { return try_nb_; }
            void set_try_nb(const unsigned int& try_nb) { try_nb_ = try_nb; }
            void increment_try_nb() { ++try_nb_; }

            /// Id
            void set_id(const size_t& id) { id_ = id; }
            size_t id() const { return id_; }

            /// Interior
            void set_is_interior(const bool& is_interior) { is_interior_ = is_interior; }
            bool is_interior() const { return is_interior_; }

            /// Operators
            bool operator==(const Self& pv) const { return ( id() == pv.id() ); }

            bool operator<(const Self& pv) const
            {
                // vertex type (smallest-interior first)
                //if ( vertex()->in_dimension() != pv.vertex()->in_dimension() )
                //    return vertex()->in_dimension() > pv.vertex()->in_dimension();
                // nb incident slivers (smallest first)
                if (is_interior() && !pv.is_interior())
                    return true;
                else if (!is_interior() && pv.is_interior())
                    return false;
                else if ( sliver_nb() != pv.sliver_nb() )
                    return sliver_nb() < pv.sliver_nb();
                // min angle (smallest first)
                else if ( min_value() != pv.min_value() )
                    return min_value() < pv.min_value();
                // try nb (smallest first)
                else if ( try_nb() != pv.try_nb() )
                    return try_nb() < pv.try_nb();
                // TODO: perturbation type (smallest first) 
                //else if ( perturbation() != pv.perturbation() )
                //    return *perturbation() < *pv.perturbation();
                return ( id() < pv.id() ); // all characteristics are the same!
            }

            /// Dummy functions
            void update_saved_erase_counter() {}
            bool is_zombie() { return false; }

            friend std::ostream& operator<<(std::ostream& os, const PVertex_& pv) {
                std::string type = pv.is_interior() ? "interior" : "boundary";
                os << "vertex("<< type <<"):" << pv.id() << "\tsliver_count:" << pv.sliver_nb() << "\tmin_value:" << pv.min_value() << std::endl;
                return os;
            }

        private:
            /// Private datas
            Vertex_handle vertex_handle_;
            unsigned int incident_sliver_nb_;
            FT min_value_;
            unsigned int try_nb_;
            PerturbationMethod p_perturbation_;
            size_t id_;
            bool is_interior_;
        };

        typedef PVertex_<FT,
            Vertex_handle,
            Point,
            PerturbationMethod> PVertex;

        class PVertex_id
        {
        public:
            typedef boost::readable_property_map_tag category;
            typedef size_t value_type;
            typedef PVertex key_type;

            value_type operator[] (const key_type& pv) const { return pv.id(); }

            friend inline
                value_type get(const PVertex_id& m, const key_type& k)
            {
                return m[k];
            }
        };

        typedef std::less<PVertex> less_PVertex;
        typedef Modifiable_priority_queue<PVertex, less_PVertex, PVertex_id> PQueue;
    
    // private function

    // public function
    public:
        Sliver_perturbation_removal(boost::shared_ptr<Triangulation> & m_tr, bool is_octree)
            : m_tr(m_tr), is_octree(is_octree)
        {}
        
        // input lower bound in degree
        int perturb(FT threshold, bool force_empty)
        {
            // found all slivers, add vertices to a priority queue
            FT tan_ub = tan(threshold / 180.0 * boost::math::constants::pi<double>());
            std::unordered_map<unsigned int, PVertex> PVertex_buffer_map;

            int count = 0;
            for (Finite_cells_iterator cb = m_tr->finite_cells_begin(); cb != m_tr->finite_cells_end(); cb++)
            {
                FT min_angle_tan = 1 / largest_cot(cb);
                if (tan_ub > min_angle_tan) // if sliver
                {
                    for (int i = 0; i < 4; i++)
                    {
                        unsigned int v_idx = cb->vertex(i)->index();
                        //if (cb->vertex(i)->type() == Triangulation::Point_type::STEINER)
                        {
                            PVertex& pv = PVertex_buffer_map[v_idx];
                            if (pv.sliver_nb() == 0)
                            {
                                pv.set_vertex(cb->vertex(i));
                                pv.set_id(v_idx);
                                pv.set_sliver_nb(1);
                                pv.set_min_value(min_angle_tan); // stores tan value corresponding to min dihedral angle
                                if (is_octree && !force_empty)
                                    pv.set_perturbation(OCTREE);
                                else
                                    pv.set_perturbation(RADIUS);
                                pv.set_is_interior(has_finite_voronoi_cell(cb->vertex(i)));
                            }
                            else
                            {
                                pv.set_sliver_nb(pv.sliver_nb() + 1);
                                if (min_angle_tan < pv.min_value())
                                    pv.set_min_value(min_angle_tan);
                            }
                        }
                    }
                    count++;
                }
                
            }
            // build queue
            PQueue pqueue(m_tr->number_of_vertices());
            for (const auto& i : PVertex_buffer_map)
                pqueue.push(i.second);
            
            // compute perturbation vector for each sliver
            while (!pqueue.empty())  
            {
                // get vector of sliver cells
                PVertex pv = pqueue.top_and_pop();
                if ((!force_empty && pv.try_nb() == 3) || (force_empty && pv.try_nb() == 5)) // avoid infinite loop
                    break; 
                int iter = 0;
                Point old_pos = pv.vertex()->point();
                unsigned char type = pv.vertex()->type();
                while (iter<10)
                {
                    iter++;
                    // compute perturbation vector
                    Sliver_info info = get_slivers(pv, tan_ub);
                    std::vector<Cell_handle> slivers = info.first;
                    if (slivers.size() == 0)
                        break;
                    Vector pertubation_vec = compute_displacement(pv, slivers, 0.1);
                    if (pertubation_vec != CGAL::NULL_VECTOR)
                    {
                        m_tr->move(pv.vertex(), old_pos + pertubation_vec);
                        // check retry condition
                        Sliver_info new_info = get_slivers(pv, tan_ub); // move method should give a same vertex handler
                        std::vector<Cell_handle> new_slivers = new_info.first;
                        if (force_empty)
                        {
                            if (new_slivers.size() != 0 || (new_info.second < info.second)) // if not 0, force slivers to be removed
                            {
                                pv.next_perturbation();
                                pv.set_vertex(m_tr->move(pv.vertex(), old_pos));
                            }
                            else
                            {
                                pv.set_min_value(new_info.second);
                                count -= (slivers.size() - new_slivers.size());
                            }
                        }
                        else 
                        {
                            //if ((new_slivers.size() > slivers.size()) || (new_info.second < info.second)) // if not better, revert and reinsert
                            if (new_info.second < info.second) // if not better, revert and reinsert
                            {   
                                pv.next_perturbation();
                                pv.set_vertex(m_tr->move(pv.vertex(), old_pos));
                            }
                            else
                            {
                                pv.set_min_value(new_info.second);
                                count -= (slivers.size() - new_slivers.size());
                            }
                        }
                    }
                    else
                        pv.next_perturbation();
                }
                if (iter == 5) // reinsert with random perturbation
                {
                    //std::vector<Cell_handle> slivers = get_slivers(pv, tan_ub).first;
                    //Vector pertubation_vec = compute_displacement(pv, slivers, 0.25); // a short random move 
                    //pv.set_vertex(m_tr->move(pv.vertex(), old_pos+pertubation_vec));
                    //std::vector<Cell_handle> new_slivers = get_slivers(pv, tan_ub).first;
                    //pv.set_sliver_nb(new_slivers.size());
                    pv.increment_try_nb();
                    pqueue.push(pv);
                }
            }
            //std::cout << "  slivers number:" << count << " with threshold " << threshold << std::endl;
            return count; 
        }

    private:
        FT cotan_per_edge(Cell_handle cell, int i, int j)
        {
            Vertex_handle vi = cell->vertex(i);
            Vertex_handle vj = cell->vertex(j);

            Point pi = vi->point();
            Point pj = vj->point();

            std::vector<Point> vpq;

            for(int i = 0; i < 4; i++)
                if(cell->vertex(i)->index() != vi->index() && cell->vertex(i)->index() != vj->index())
                    vpq.push_back(cell->vertex(i)->point());

            Vector ni = CGAL::cross_product(pi - vpq[0], pi - vpq[1]);
            Vector nj = CGAL::cross_product(pj - vpq[0], pj - vpq[1]);

            ni = ni / std::sqrt(ni * ni);
            nj = nj / std::sqrt(nj * nj);

            Vector nij = CGAL::cross_product(ni, nj);
            FT cotan = (ni * nj) / std::sqrt(nij * nij);

            return cotan;
        }

        FT largest_cot(Cell_handle cell)
        {
            FT max_cotan = -1e7;
            for(int i = 0; i < 3; i++)
                for (int j = i + 1; j < 4; j++)
                {
                    double cotan = cotan_per_edge(cell, i, j);
                    if(cotan > max_cotan) max_cotan = cotan;
                }
            return max_cotan;
        };

        bool has_finite_voronoi_cell(Vertex_handle v)
        {
            std::list<Cell_handle> cells;
            m_tr->incident_cells(v, std::back_inserter(cells));

            if(cells.size() == 0)
                return false;

            typename std::list<Cell_handle>::iterator it;
            for(it = cells.begin(); it != cells.end(); it++)
            {
                if(m_tr->is_infinite(*it))
                    return false;
            }
            return true;
        }

        int update_priority_queue(const PVertex& pv, PQueue& pqueue) const
        {
            if ( pqueue.contains(pv) )
            {
                if ( pv.is_perturbable() )
                {
                    pqueue.update(pv);
                    return 0;
                }
                else
                {
                    pqueue.erase(pv);
                    return -1;
                }
            }
            else
            {
                if ( pv.is_perturbable() )
                {
                    pqueue.push(pv);
                    return 1;
                }
            }

            return 0;
        }

        FT edge_sq_length(const Edge& e)
        {
            typename Geom_traits::Compute_squared_distance_3 sq_distance =
                m_tr->geom_traits().compute_squared_distance_3_object();

            const Point& p = m_tr->point(e.first, e.second);
            const Point& q = m_tr->point(e.first, e.third);

            return sq_distance(p,q);
        }

        FT min_incident_edge_sq_length(const Vertex_handle& v)
        {
            CGAL_precondition(!m_tr->is_infinite(v));

            // Get all incident edges
            std::vector<Edge> edges;
            m_tr->finite_incident_edges(v, std::back_inserter(edges));
            CGAL_assertion(!edges.empty());

            // Get squared min length
            typename std::vector<Edge>::iterator eit = edges.begin();
            FT min_sq_length = edge_sq_length(*eit++);

            for ( ; eit != edges.end() ; ++eit )
            {
                min_sq_length = (std::min)(min_sq_length, edge_sq_length(*eit));
            }

            return min_sq_length;
        }

        Sliver_info get_slivers(const PVertex& pv, FT tan_ub)
        {
            FT min_of_min = 1e7;
            std::vector<Cell_handle> cells, slivers;
            m_tr->finite_incident_cells(pv.vertex(), std::back_inserter(cells));
            for (auto c : cells) { // extract slivers into a vector
                double min_angle_tan = 1 / largest_cot(c);
                if (tan_ub > min_angle_tan)
                { 
                    if (min_of_min > min_angle_tan)
                        min_of_min = min_angle_tan;
                    slivers.push_back(c);
                }
            }
            return Sliver_info(slivers, min_of_min);
        }
        
        Vector compute_displacement(const PVertex& pv, std::vector<Cell_handle> slivers, FT alpha) 
        {
            Vector grad;
            if (slivers.size() == 0)
                return CGAL::NULL_VECTOR;
            else if (slivers.size() > 2)
            {
                alpha = 0.25;
                grad = compute_random_perturbation();
            }
            else if (pv.perturbation() == OCTREE)
            {   
                grad = compute_octree_perturbation(slivers);
                alpha = 0.25;
            }
            else if (pv.perturbation() == RADIUS)
                grad = compute_gradient_radius(pv.vertex(), slivers);
            else if (pv.perturbation() == VOLUME)
                grad = compute_gradient_volume(pv.vertex(), slivers);
            else if (pv.perturbation() == RANDOM)
            { 
                alpha = 0.25;
                grad = compute_random_perturbation();
            }
            FT grad_length = sqrt(grad * grad);
            FT norm_length = sqrt(min_incident_edge_sq_length(pv.vertex()));
            grad = grad / grad_length * norm_length * alpha;
            if (isnan(grad.x()) || isnan(grad.y()) || isnan(grad.z()))
                return CGAL::NULL_VECTOR;
            else
                return grad;
        }

        Vector compute_octree_perturbation(const std::vector<Cell_handle>& slivers) const
        {
            switch (slivers.size())
            {
            case 1:
                return compute_octree_perturbation_(slivers.front());
                break;
            case 2:
                {
                    Vector v1 = compute_octree_perturbation_(slivers.front());
                    Vector v2 = compute_octree_perturbation_(slivers.back());
                    if (v1 * v2 > 0)
                        return 0.5 * (v1 + v2);
                    else
                        return 0.5 * (v1 - v2); 
                    break;
                }
            default:
                break;
            }

            // May happen if sq_radius_gradient is not relevant for this vertex
            return CGAL::NULL_VECTOR;
        }

        // for octree plain sliver, we know where to go (at least somehow)
        Vector compute_octree_perturbation_(const Cell_handle& cell) const 
        {
            Point v[3];
            for (int i = 0; i < 3; i++)
            {
                v[i] = cell->vertex(i)->point();
            }
            Vector pert = CGAL::cross_product(v[2] - v[0], v[1] - v[0]); 
            return pert;
        }
        
        Vector compute_random_perturbation()
        {
            typedef boost::lagged_fibonacci607 base_generator_type;
            static base_generator_type generator_;
            static boost::uniform_real<FT> uni_dist_;
            static boost::variate_generator<base_generator_type&,
                boost::uniform_real<FT> > random_(generator_, uni_dist_);
            Vector rnd_vector(random_() - 0.5, random_()-0.5, random_()-0.5);
            return rnd_vector;
        }

        Vector compute_gradient_radius(
            const Vertex_handle& v,
            const std::vector<Cell_handle>& slivers) const
        {
            switch (slivers.size())
            {
            case 1:
                return compute_gradient_radius_(v, slivers.front());
                break;
            case 2:
            {
                Vector v1 = compute_gradient_radius_(v, slivers.front());
                Vector v2 = compute_gradient_radius_(v, slivers.back());
                if( v1 * v2 > 0 )
                    // "+0.5" because sq_radius has to go up
                    return 0.5*(v1 + v2);
                break;
            }
            default:
                break;
            }

            // May happen if sq_radius_gradient is not relevant for this vertex
            return CGAL::NULL_VECTOR;
        }

        Vector compute_gradient_radius_(const Vertex_handle& v, const Cell_handle& cell) const
        {
            typename Gt::Construct_translated_point_3 translate =
                m_tr->geom_traits().construct_translated_point_3_object();

            unsigned int index = cell->index(v);

            const Point& wvp = m_tr->point(cell, index);
            const Point& wp2 = m_tr->point(cell, (index+1)&3);
            const Point& wp3 = m_tr->point(cell, (index+2)&3);
            const Point& wp4 = m_tr->point(cell, (index+3)&3);

            // translate the tet so that 'wp4' is the origin
            Vector translate_to_origin(CGAL::ORIGIN, wp4);
            const Point& p1 = translate(wvp, - translate_to_origin);
            const Point& p2 = translate(wp2, - translate_to_origin);
            const Point& p3 = translate(wp3, - translate_to_origin);

            // pre-compute everything
            FT sq_p1 = p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z();
            FT sq_p2 = p2.x()*p2.x() + p2.y()*p2.y() + p2.z()*p2.z();
            FT sq_p3 = p3.x()*p3.x() + p3.y()*p3.y() + p3.z()*p3.z();

            // every derivative is computed w.r.t p1 (x1, y1, z1)
            FT da_dx = p2.y()*p3.z() - p3.y()*p2.z();
            FT da_dy = p2.z()*p3.x() - p2.x()*p3.z();
            FT da_dz = p2.x()*p3.y() - p3.x()*p2.y();

            FT dDx_dx = -2*p1.x()*da_dx;
            FT dDx_dy = -2*p1.y()*da_dx + sq_p2*p3.z() - sq_p3*p2.z();
            FT dDx_dz = -2*p1.z()*da_dx - sq_p2*p3.y() + sq_p3*p2.y();

            FT dDy_dx = -2*p1.x()*da_dy - sq_p2*p3.z() + sq_p3*p2.z();
            FT dDy_dy = -2*p1.y()*da_dy;
            FT dDy_dz = -2*p1.z()*da_dy + sq_p2*p3.x() - sq_p3*p2.x();

            FT dDz_dx = -2*p1.x()*da_dz + sq_p2*p3.y() - sq_p3*p2.y();
            FT dDz_dy = -2*p1.y()*da_dz - sq_p2*p3.x() + sq_p3*p2.x();
            FT dDz_dz = -2*p1.z()*da_dz;

            FT a  = p1.x()*da_dx + p1.y()*da_dy + p1.z()*da_dz;
            if ( CGAL_NTS is_zero(a) )
                return CGAL::NULL_VECTOR;

            FT Dx = -sq_p1*da_dx + p1.y()*(sq_p2*p3.z() - sq_p3*p2.z()) - p1.z()*(sq_p2*p3.y() - sq_p3*p2.y());
            FT Dy = -sq_p1*da_dy - p1.x()*(sq_p2*p3.z() - sq_p3*p2.z()) + p1.z()*(sq_p2*p3.x() - sq_p3*p2.x());
            FT Dz = -sq_p1*da_dz + p1.x()*(sq_p2*p3.y() - sq_p3*p2.y()) - p1.y()*(sq_p2*p3.x() - sq_p3*p2.x());

            // compute gradient vector
            FT sum_sqD = Dx*Dx + Dy*Dy + Dz*Dz;
            FT gx = (Dx*dDx_dx + Dy*dDy_dx + Dz*dDz_dx) / (2.0*a*a) - (da_dx * sum_sqD) / (2.0*a*a*a);
            FT gy = (Dx*dDx_dy + Dy*dDy_dy + Dz*dDz_dy) / (2.0*a*a) - (da_dy * sum_sqD) / (2.0*a*a*a);
            FT gz = (Dx*dDx_dz + Dy*dDy_dz + Dz*dDz_dz) / (2.0*a*a) - (da_dz * sum_sqD) / (2.0*a*a*a);

            return Vector(gx, gy, gz);
        };

        Vector compute_gradient_volume( const Vertex_handle& v, const std::vector<Cell_handle>& slivers) const
        {
            switch (slivers.size())
            {
            case 1:
                return -1*compute_gradient_volume_(v, slivers.front());
                break;
            case 2:
            {
                Vector v1 = compute_gradient_volume_(v, slivers.front());
                Vector v2 = compute_gradient_volume_(v, slivers.back());
                if( v1 * v2 > 0 )
                    // "-0.5" because volume has to go down
                    return -0.5 * (v1 + v2);
                break;
            }
            default:
                break;
            }

            // May happen if volume_gradient is not relevant for this vertex
            return CGAL::NULL_VECTOR;
        }

        Vector compute_gradient_volume_(const Vertex_handle& v, const Cell_handle& cell) const
        {
            CGAL_precondition(cell->has_vertex(v));

            const int i = cell->index(v);

            // fixed vertices: (the ones with index != i)
            int k1 = (i+1)&3;
            int k2 = (i+2)&3;
            int k3 = (i+3)&3;

            if ( (i&1) == 0 )
                std::swap(k1,k3);

            const Point& p1 = m_tr->point(cell, k1);
            const Point& p2 = m_tr->point(cell, k2);
            const Point& p3 = m_tr->point(cell, k3);

            FT gx =  p2.y()*p3.z() + p1.y()*(p2.z()-p3.z())
                - p3.y()*p2.z() - p1.z()*(p2.y()-p3.y());

            FT gy = -p2.x()*p3.z() - p1.x()*(p2.z()-p3.z())
                + p3.x()*p2.z() + p1.z()*(p2.x()-p3.x());

            FT gz =  p2.x()*p3.y() + p1.x()*(p2.y()-p3.y())
                - p3.x()*p2.y() - p1.y()*(p2.x()-p3.x());

            return (1.0 / 6.0 * Vector(gx, gy, gz));
        }
        
    };


    // wrapper
    template<class Gt,
        class Tr>
    bool remove_sliver(boost::shared_ptr<Tr>& tr, double threshold, bool is_octree)
    {
        Sliver_perturbation_removal<Gt, Tr> perturber(tr, is_octree);
        if (perturber.perturb(threshold, false) > 0)
            if (perturber.perturb(threshold, true) > 0)
                perturber.perturb(1, true); // to get rid of exact plain slivers that are very harmful
        return true;
    }

}
#endif // CGAL_IMPLICIT_RECONSTRUCTION_PERTURB_SLIVER_H