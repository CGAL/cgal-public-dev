
namespace CGAL {

enum Sides {
    ARR_LEFT, ARR_RIGHT, ARR_TOP, ARR_BOTTOM
};

class Boundary_type_tag {};

class Bounded_tag : public Boundary_type_tag {};
class Unbounded_tag : public Boundary_type_tag {};
class Contraction_tag : public Boundary_type_tag {};
class Identification_tag : public Boundary_type_tag {};

enum Boundary {
    BOUNDED,
    UNBOUNDED_FICT,
    UNBOUNDED_VERTEX,
    CONTRACTION,
    IDENTIFICATION
};


#define CGAL_ARR_SIDE_TYPEDEFS \
    typedef typename Dcel::Vertex Vertex; \
    typedef typename Dcel::Halfedge Halfedge; \
    typedef typename Dcel::Face Face; \
    \
    typedef typename Geo_traits_2::Point_2 Point_2; \
    typedef typename Geo_traits_2::X_monotone_curve_2 X_monotone_curve_2; \

// end


template < class Dcel_, class GeoTraits_2 >
class Side {
public:
    typedef Dcel_ Dcel;
    typedef GeoTraits_2 Geo_traits_2;

    CGAL_ARR_SIDE_TYPEDEFS;

public:
    Side(Dcel *dcel) :
        _m_dcel(dcel) {
    }

    virtual ~Side() {
        delete _m_dcel;
    }
    
    virtual void update(const Point_2& point) {
        std::cout << "TT-Side: update point" << std::endl;
    }

    virtual void update(const X_monotone_curve_2& xc, 
                        const CGAL::Curve_end& ind) {
        std::cout << "TT-Side: update cve" << std::endl;
    }
    
protected:
    Dcel *_m_dcel;
};

template < class Dcel_, class GeoTraits_2 >
class Fictitious_halfedge_side : public Side< Dcel_, GeoTraits_2 > {

public:
    typedef Dcel_ Dcel;
    typedef GeoTraits_2 Geo_traits_2;

    typedef Side< Dcel, Geo_traits_2 > Base;

    CGAL_ARR_SIDE_TYPEDEFS;
    
    Fictitious_halfedge_side(Dcel *dcel) :
        Base(dcel) {
        construct();
    }
    
    virtual ~Fictitious_halfedge_side() {
        delete _m_min_vertex;
        delete _m_max_vertex;
    }
    
    virtual void update(const Point_2& point) {
        std::cout << "TT-FH-Side: update point" << std::endl;
    }

    virtual void update(const X_monotone_curve_2& xc, 
                        const CGAL::Curve_end& ind) {
        std::cout << "TT-FH-Side: update cve" << std::endl;
    }
    
private:
    
    void construct() {
        _m_min_vertex = this->_m_dcel->new_vertex();
        _m_max_vertex = this->_m_dcel->new_vertex();
        
        // side v_bl->set_boundary (MINUS_INFINITY, MINUS_INFINITY);
          
        Halfedge *he = this->_m_dcel->new_edge();
        Halfedge *he_t = he->opposite();
        
        he->set_curve (NULL);
        he_t->set_curve (NULL);

        // Set the direction of the halfedges:
        he->set_direction (LEFT_TO_RIGHT); // TODO

        he->set_vertex (_m_min_vertex);     he_t->set_vertex (_m_max_vertex);

        _m_min_vertex->set_halfedge (he);
        _m_max_vertex->set_halfedge (he_t);
    }
    
protected:
    Vertex *_m_min_vertex;
    Vertex *_m_max_vertex;
    
}; 

template < class Dcel_, class GeoTraits_2 >
class Vertex_side : public Side< Dcel_, GeoTraits_2 > {
    
public:
    typedef Dcel_ Dcel;
    typedef GeoTraits_2 Geo_traits_2;

    typedef Side< Dcel, Geo_traits_2 > Base;

    CGAL_ARR_SIDE_TYPEDEFS;

    Vertex_side(Dcel *dcel) : 
        Base(dcel),
        _m_vertex(NULL) {
    }

    virtual ~Vertex_side() {
        delete _m_vertex;
    }
    
    virtual void update(const Point_2& point) {
        std::cout << "TT-V-Side: update point" << std::endl;
    }

    virtual void update(const X_monotone_curve_2& xc, 
                        const CGAL::Curve_end& ind) {
        std::cout << "TT-V-Side: update cve" << std::endl;
    }

protected:
    Vertex *_m_vertex;
};

template < class Dcel_, class GeoTraits_2 >
class Identification : public Side< Dcel_, GeoTraits_2 > {
    
public:
    typedef Dcel_ Dcel;
    typedef GeoTraits_2 Geo_traits_2;

    typedef Side< Dcel, Geo_traits_2 > Base;

    CGAL_ARR_SIDE_TYPEDEFS;

    Identification(Dcel *dcel) : 
        Base(dcel) {
    }
    
    virtual ~Identification() {

    }
    

    virtual void update() {
    }

    
};


template < class Dcel_, class GeoTraits_2 >
class Topology_traits {
public:

    typedef Dcel_ Dcel;
    typedef GeoTraits_2 Geo_traits_2;
    
    typedef Side< Dcel, Geo_traits_2 > Side;
    typedef Vertex_side< Dcel, Geo_traits_2 > Vertex_side;
    typedef Fictitious_halfedge_side< Dcel, Geo_traits_2 > 
    Fictitious_halfedge_side;
    typedef Identification< Dcel, Geo_traits_2 > Identification;    

    Topology_traits(
            Boundary left,
            Boundary bottom,
            Boundary top,
            Boundary right
    ) {

        CGAL_precondition(
                left != IDENTICATION && right != IDENTIFICATION ||
                left == IDENTICATION && right == IDENTIFICATION
        );

        CGAL_precondition(
                top != IDENTICATION && bottom != IDENTIFICATION ||
                top == IDENTICATION && bottom == IDENTIFICATION
        );

        // TODO initial interior "face" + boundary
#if 0
        Outer_ccb          *oc = this->m_dcel.new_outer_ccb();
        Inner_ccb          *ic = this->m_dcel.new_inner_ccb();
        Face               *in_f = this->m_dcel.new_face();

        this->m_dcel.new_face();
        
        fict_face->set_unbounded (true);
        fict_face->set_fictitious (true);

        oc->set_face (in_f);
        ic->set_face (fict_face);
        
        he1->set_inner_ccb (ic);       he1_t->set_outer_ccb (oc);
        he2->set_inner_ccb (ic);       he2_t->set_outer_ccb (oc);
        he3->set_inner_ccb (ic);       he3_t->set_outer_ccb (oc);
        he4->set_inner_ccb (ic);       he4_t->set_outer_ccb (oc);

        // Set the inner component of the fictitious face.
        fict_face->add_inner_ccb (ic, he1);
        
        // Set the real unbounded face, in the interior of the bounding rectangle.
        in_f->add_outer_ccb (oc, he1_t);
        in_f->set_unbounded (true);

        n_inf_vertex = 4;
#endif
        
        switch (left) {
        case UNBOUNDED_FICT:
            this->_m_left = new Fictitious_halfedge_side(&this->_m_dcel);
            break;
        case UNBOUNDED_VERTEX:
        case CONTRACTION:
            this->_m_left = new Vertex_side(&this->_m_dcel);
            break;
        case IDENTIFICATION:
            this->_m_left = new Identification(&this->_m_dcel);
        default:
            break;
        }
        
        switch (bottom) {
        case UNBOUNDED_FICT:
            this->_m_bottom = new Fictitious_halfedge_side(&this->_m_dcel);
            break;
        case UNBOUNDED_VERTEX:
        case CONTRACTION:
            this->_m_bottom = new Vertex_side(&this->_m_dcel);
            break;
        case IDENTIFICATION:
            this->_m_bottom = new Identification(&this->_m_dcel);
        default:
            break;
        }
        
        switch (top) {
        case UNBOUNDED_FICT:
            this->_m_top = new Fictitious_halfedge_side(&this->_m_dcel);
            break;
        case UNBOUNDED_VERTEX:
        case CONTRACTION:
            this->_m_top = new Vertex_side(&this->_m_dcel);
            break;
        case IDENTIFICATION:
            this->_m_top = new Identification(&this->_m_dcel);
        default:
            break;
        }
        
        switch (right) {
        case UNBOUNDED_FICT:
            this->_m_right = new Fictitious_halfedge_side(&this->_m_dcel);
            break;
        case UNBOUNDED_VERTEX:
        case CONTRACTION:
            this->_m_right = new Vertex_side(&this->_m_dcel);
            break;
        case IDENTIFICATION:
            this->_m_right = new Identification(&this->_m_dcel);
        default:
            break;
        }

        // TODO connect the sides
        
    }
    
public:
    
        
    void query() {} // locate curve-end, pt,

    
    // TODO boundary_type can be answered after construction


private:
    Dcel _m_dcel;

    Side *_m_left;
    Side *_m_bottom;
    Side *_m_top;
    Side *_m_right;

};


} // namespace CGAL
