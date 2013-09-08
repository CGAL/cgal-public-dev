
#ifndef files_Arr_cmap_dcel_h
#define files_Arr_cmap_dcel_h
/*! \file
 * The definition of the Arr_cmap_dcel<Traits, CMap> class.
 */

#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/Combinatorial_map_with_holes.h>
#include <list>

namespace CGAL {
/*! \class
 * The arrangement DCEL class represented by Combinatorial_map_with_holes.
 * The Traits parameters corresponds to a geometric traits class, which
 * defines the Point_2 and X_monotone_curve_2 types.
 */

//Forward declaration
template <class Traits_, class CMap> class Arr_cmap_vertex;
template <class Traits_, class CMap> class Arr_cmap_halfedge;
template <class Traits_, class CMap> class Arr_cmap_face;
template <class Traits_, class CMap> class Arr_cmap_outer_ccb;
template <class Traits_, class CMap> class Arr_cmap_inner_ccb;
template <class Traits_, class CMap> class Arr_cmap_isolated_vertex;
template <class Traits_, class Cmap> class Arr_cmap_dcel;

/*! \class
 * The arrangement DCEL vertex class.
 */
template <class Traits_, class CMap>
class Arr_cmap_vertex
{
public:

  typedef Traits_                                   Traits;
  typedef typename Traits::Point_2                  Point;
  typedef typename CMap::Dart_handle                Dart_handle;

public:

  typedef Arr_cmap_vertex<Traits, CMap>             Vertex;
  typedef Arr_cmap_halfedge<Traits, CMap>           Halfedge;
  //typedef Arr_cmap_face<Traits, CMap>               Face;
  typedef Arr_cmap_isolated_vertex<Traits, CMap>    Isolated_vertex;
  typedef Arr_cmap_dcel<Traits, CMap>               Dcel;

public:

  //Dcel* dcel;
  Dart_handle dh;

private:

  Point           *pt;     //The point associated with the vertex.
  Halfedge        *inc;    //The incident halfedge.
  Isolated_vertex *civ;    //The isolated vertex information.
  char            pss[2];  //The x and y parameter spaces
                           //(condensed in two bytes).

public:

  //default constructor
  Arr_cmap_vertex () :
    pt(NULL),
    inc(NULL),
    civ(NULL)
    //dcel(NULL)
  {
    pss[0] = pss[1] = static_cast<char> (CGAL::ARR_INTERIOR);
    dh = CMap::null_dart_handle;
  }

  /*! Destructor. */
  virtual ~Arr_cmap_vertex () {}

  /*! Check if the vertex is isolated. */
  bool is_isolated () const
  {
    return (dh == CMap::null_dart_handle);
  }

  /*! Get an incident halfedge (const version). */
  const Halfedge* halfedge() const
  {
    CGAL_precondition(! is_isolated());
    return this->inc;
  }

  /*! Get an incident halfedge (non-const version). */
  Halfedge* halfedge()
  {
    CGAL_precondition(! is_isolated());
    return this->inc;
  }

  /*! Get the isolated vertex information (non-const version). */
  Isolated_vertex* isolated_vertex()
  {
    CGAL_precondition(is_isolated());
    return this->civ;
  }

  /*! Get the isolated vertex information (const version). */
  const Isolated_vertex* isolated_vertex() const
  {
    CGAL_precondition(is_isolated());
    return this->civ;
  }

  /*! Set the isolated vertex information. */
  void set_isolated_vertex(Isolated_vertex* iv)
  {
    civ = iv;
  }

  /*! Get the point (const version). */
  const Point& point() const
  {
    CGAL_assertion (pt != NULL);
    return (*pt);
  }

  /*! Get the point (non-const version). */
  Point& point()
  {
    CGAL_assertion (pt != NULL);
    return (*pt);
  }

  /*! Set the point (may be a NULL point). */
  void set_point (Point *p)
  {
    pt = p;
  }

  /*! Check if the point pointer is NULL. */
  bool has_null_point () const
  {
    return (pt == NULL);
  }

  /*! Set an incident halfedge (for non-isolated vertices). */
  void set_halfedge(Halfedge* e)
  {
    this->inc = e;
    dh = e->dh;
  }

  /*! Get the boundary type in x. */
  Arr_parameter_space parameter_space_in_x() const
  {
    return (Arr_parameter_space (pss[0]));
  }

  /*! Get the boundary type in y. */
  Arr_parameter_space parameter_space_in_y() const
  {
    return (Arr_parameter_space (pss[1]));
  }

  /*! Set the boundary conditions of the vertex. */
  void set_boundary (Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    pss[0] = static_cast<char> (ps_x);
    pss[1] = static_cast<char> (ps_y);
    return;
  }

  /*! Assign from another vertex. */
  virtual void assign (const Arr_cmap_vertex<Traits, CMap>& v)
  {
    pt = v.pt;
    pss[0] = v.pss[0];
    pss[1] = v.pss[1];
  }
};

/*! \class
 * The arrangement DCEL halfedge class.
 */
template <class Traits_, class CMap>
class Arr_cmap_halfedge
{
public:

  typedef Traits_                                    Traits;
  typedef typename Traits::X_monotone_curve_2        X_monotone_curve;
  typedef typename CMap::Dart_handle                 Dart_handle;

public:

  typedef Arr_cmap_vertex<Traits, CMap>              Vertex;
  typedef Arr_cmap_halfedge<Traits, CMap>            Halfedge;
  typedef Arr_cmap_face<Traits, CMap>                Face;
  typedef Arr_cmap_hole<CMap>                        Inner_ccb;
  typedef Arr_cmap_outer_ccb<Traits, CMap>           Outer_ccb;
  typedef Arr_cmap_dcel<Traits, CMap>                Dcel;

public:

  Dcel* dcel;         // The corresponding DCEL class
  Dart_handle dh;     // The correspoinding Dart_handle

private:

  X_monotone_curve *xv;  // The associated x-monotone curve.
  Halfedge *oppo;        // The opposite halfedge.
  Halfedge *pre;         // The previous halfedge in the component boundary.
  Halfedge *nex;         // The next halfedge in the component boundary.
  Vertex *v;             // The incident vertex (the target of the halfedge).
  bool left_to_right;    // Determine the direction of the halfedge
  bool on_inner_ccb;     // Determine whether it is on inner ccb
  Outer_ccb  *outer;
  Inner_ccb *inner;
    
public:

  /*! Default constructor */
  Arr_cmap_halfedge() :
    xv(NULL),
    v(NULL),
    oppo(NULL),
    pre(NULL),
    nex(NULL),
    dcel(NULL)
  {
    dh = CMap::null_dart_handle;
  }

  /*! Destructor. */
  virtual ~Arr_cmap_halfedge()
  {}

  /*! Get the opposite halfedge (non-const version). */
  Halfedge* opposite()
  {
    return this->oppo;
  }

  /*! Get the opposite halfedge (const version). */
  const Halfedge* opposite() const
  {
    return this->oppo;
  }

  /*! Get the previous halfedge along the chain (non-const version). */
  Halfedge* prev()
  {
    return this->pre;
  }

  /*! Get the previous halfedge along the chain (const version). */
  const Halfedge* prev() const
  {
    return this->pre;
  }

  /*! Get the next halfedge along the chain (non-const version). */
  Halfedge* next()
  {
    return this->nex;
  }

  /*! Get the next halfedge along the chain (const version). */
  const Halfedge* next() const
  {
    return this->nex;
  }

  /*! Get the x-monotone curve (non-const version). */
  X_monotone_curve& curve()
  {
    CGAL_precondition(xv != NULL);
    return (*xv);
  }

  /*! Get the x-monotone curve (const version). */
  const X_monotone_curve& curve() const
  {
    CGAL_precondition(xv != NULL);
    return (*xv);
  }

  /*! Check if the curve pointer is NULL. */
  bool has_null_curve() const
  {
    return (xv == NULL);
  }

  /*! Sets the opposite halfedge. */
  void set_opposite(Halfedge* opp)
  {
    this->oppo = opp;
  }

  /*! Get the direction of the halfedge. */
  Arr_halfedge_direction direction () const
  {
    if (left_to_right)
      return (ARR_LEFT_TO_RIGHT);
    else
      return (ARR_RIGHT_TO_LEFT);
  }

  /*! Set the direction of the edge (and of its opposite halfedge). */
  void set_direction() const
  {
    if (dir == ARR_LEFT_TO_RIGHT)
    {
      left_to_right = true;
      oppo.left_to_right = false;
    }
    else
    {
      left_to_right = false;
      oppo.left_to_right = true;
    }
  }

  /*! Set the previous halfedge along the chain. */
  void set_prev(Halfedge* prev)
  {
    this->pre = prev;
    prev->nex = this;
    dcel->cm.template link_beta<1>(prev->dh, dh);
  }

  /*! Set the next halfedge along the chain. */
  void set_next(Halfedge* next)
  {
    this->nex = next;
    next->pre = this;
    dcel->cm.template link_beta<1>(dh, nex->dh);
  }

  /*! Set the target vertex. */
  void set_vertex(Vertex* v)
  {
    this->v = v;
    v->dh = dh;
  }

  /*! Set the x-monotone curve. */
  void set_curve(X_monotone_curve* c)
  {
    xv = c;
    oppo->xv = c;
  }

  /*! Assign from another halfedge. */
  virtual void assign(const Arr_cmap_halfedge<Traits, CMap>& he)
  {
    xv = he.xv;
  }

  /*! Get the target vertex (const version). */
  const Vertex* vertex () const
  {
    return (this->v);
  }

  /*! Get the target vertex (non-const version). */
  Vertex* vertex () const
  {
    return (this->v);
  }

  /*! Check whether the halfedge lies on the boundary of an inner CCB. */
  bool is_on_inner_ccb () const
  {
    return (on_inner_ccb);
  }

  /*!
   * Get an incident outer CCB (const version).
   * \pre The edge does not lie on an inner CCB.
   */
  const Outer_ccb* outer_ccb () const
  {
    CGAL_precondition (! is_on_inner_ccb());
    return this->outer;
  }

  /*!
   * Get an incident outer CCB (non-const version).
   * \pre The edge does not lie on an inner CCB.
   */
  Outer_ccb* outer_ccb ()
  {
    CGAL_precondition (! is_on_inner_ccb());
    return this->outer;
  }

  /*! Set the incident outer CCB. */
  void set_outer_ccb (Outer_ccb *oc)
  {
    this->outer = oc;
    on_inner_ccb = false;
  }

  /*!
   * Get an incident inner CCB (const version).
   * \pre The edge lies on an inner CCB.
   */
  const Inner_ccb* inner_ccb () const
  {
    CGAL_precondition (is_on_inner_ccb());
    return this->inner;
  }

  /*!
   * Get an incident inner CCB (non-const version).
   * \pre The edge lies on an inner CCB.
   */
  Inner_ccb* inner_ccb ()
  {
    CGAL_precondition (is_on_inner_ccb());
    return this->inner;
  }

  /*! Set the incident inner CCB. */
  void set_inner_ccb (Inner_ccb *ic)
  {
    this->inner = ic;
    on_inner_ccb = true;
  }
};

/*! \class
 * The arrangement DCEL face class.
 */
template <class Traits_, class CMap>
class Arr_cmap_face
{
public:

  typedef Traits_                        Traits;
  typedef typename CMap::Dart_handle     Dart_handle;

public:

  typedef Arr_cmap_vertex<Traits, CMap>              Vertex;
  typedef Arr_cmap_halfedge<Traits, CMap>            Halfedge;
  typedef Arr_cmap_face<Traits, CMap>                Face;
  typedef Arr_cmap_inner_ccb<Traits, CMap>           Inner_ccb;
  typedef Arr_cmap_outer_ccb<Traits, CMap>           Outer_ccb;
  typedef Arr_cmap_isolated_vertex<Traits, CMap>     Isolated_vertex;
  typedef Arr_cmap_dcel<Traits, CMap>                Dcel;

public:

  typedef std::list<Outer_ccb>                Outer_ccbs_container;
  typedef Outer_ccbs_container::iterator      Outer_ccb_iterator;
  typedef Outer_ccbs_container::const_iterator
                                              Outer_ccb_const_iterator;

  typedef std::list<Inner_ccb>                Inner_ccbs_container;
  typedef Inner_ccbs_container::iterator      Inner_ccb_iterator;
  typedef Inner_ccbs_container::const_iterator
                                              Inner_ccb_const_iterator;

  typedef std::list<Isolated_vertex>          Isolated_vertex_container;
  typedef Isolated_vertex_container::iterator
                                              Isolated_vertex_iterator;
  typedef Isolated_vertex_container::const_iterator
                                              Isolated_vertex_const_iterator;

public:

  typedef CMap::Face           F;
  
  Dcel                         *dcel;
  F                            f;
  Inner_ccbs_container         inner_ccbs;
  Outer_ccbs_container         outer_ccbs;
  Isolated_vertices_container  iso_verts;

    int flag;
public:

  /*! Default constructor. */
  Arr_cmap_face() :
    dcel(NULL),
    f(NULL)
  {
    flag = 0;
  }
    
  /*! Destructor. */
  virtual ~Arr_cmap_face()
  {}

  bool is_unbounded()
  {
    return (f.dart_list.front() == CMap::null_dart_handle);
  }

  void set_unbounded(bool flag)
  {
    
  }





  /*! Get the number of outer CCBs the face has. */
  size_t number_of_outer_ccbs () const
  {
    return outer_ccbs.size();
  }

  /*! Get an iterator for the first outer CCB of the face. */
  Outer_ccb_iterator outer_ccbs_begin()
  {
    return outer_ccbs.begin();
  }

  /*! Get a past-the-end iterator for the outer CCBs inside the face. */
  Outer_ccb_iterator outer_ccbs_end()
  {
    return outer_ccbs.end();
  }

  /*! Get an const iterator for the first outer CCB inside the face. */
  Outer_ccb_const_iterator outer_ccbs_begin() const
  {
    return outer_ccbs.begin();
  }

  /*! Get a const past-the-end iterator for the outer CCBs inside the face. */
  Outer_ccb_const_iterator outer_ccbs_end() const
  {
    return outer_ccbs.end();
  }

  /*! Add an outer CCB to the face. */
  void add_outer_ccb (Outer_ccb *oc, Halfedge *h)
  {
    oc->set_iterator(h);
    f.add_boundary(h->dh);
    outer_ccbs.push_back(*oc);
  }

  /*! Erase an outer CCB of the face. */
  void erase_outer_ccb (Outer_ccb *oc)
  {
    outer_ccbs.remove(*oc);
    f.erase_boundary();
    //oc->dh = CMap::null_dart_handle;
  }








  /*! Get the number of inner CCBs the face has. */
  size_t number_of_inner_ccbs() const
  {
    return Inner_ccbs.size();
  }

  /*! Get an iterator for the first inner CCB of the face. */
  Inner_ccb_iterator inner_ccbs_begin()
  {
    return Inner_ccbs.begin();
  }

  /*! Get an const iterator for the first inner CCB inside the face. */
  Inner_ccb_const_iterator inner_ccbs_begin() const
  {
    return Inner_ccbs.begin();
  }

  /*! Get a past-the-end iterator for the inner CCBs inside the face. */
  Inner_ccb_iterator inner_ccbs_end()
  {
    return Inner_ccbs.end();
  }

  /*! Get a const past-the-end iterator for the inner CCBs inside the face. */
  Inner_ccb_const_iterator inner_ccbs_end() const
  {
    return Inner_ccbs.end();
  }

  /*! Add an inner CCB to the face. */
  void add_inner_ccb (Inner_ccb *ic, Halfedge *h)
  {
    ic->set_iterator(h);
    f.add_hole(h->dh);
    inner_ccbs.push_back(*ic);
  }

  /*! Erase an inner CCB of the face. */
  void erase_inner_ccb (Inner_ccb *ic)
  {
    inner_ccbs.remove(*ic);
    f.erase_hole(ic->dh);
  }

  // Backward compatibility:
  size_t number_of_holes () const { return number_of_inner_ccbs(); }
  Hole_iterator holes_begin() { return inner_ccbs_begin(); }
  Hole_iterator holes_end() { return inner_ccbs_end(); }
  Hole_const_iterator holes_begin() const { return inner_ccbs_begin(); }
  Hole_const_iterator holes_end() const { return inner_ccbs_end(); }





  /*! Get the number of isloated vertices inside the face. */
  size_t number_of_isolated_vertices() const
  {
    return Isolated_vertices.size();
  }

  /*! Get an iterator for the first isloated vertex inside the face. */
  Isolated_vertex_iterator isolated_vertices_begin()
  {
    return isolated_vertices.begin();
  }

  /*! Get an const iterator for the first isloated vertex inside the face. */
  Isolated_const_vertex_iterator isolated_vertices_begin() const
  {
    return isolated_vertices.begin();
  }

  /*! Get a past-the-end iterator for the isloated vertices inside the face. */
  Isolated_vertex_iterator isolated_vertices_end()
  {
    return isolated_vertices.end();
  }

  /*! Get a const past-the-end iterator for the isloated vertices inside the
   * face. */
  Isolated_const_vertex_iterator isolated_vertices_end() const
  {
    return isolated_vertices.end();
  }

  /*! Add an isloated vertex inside the face. */
  void add_isolated_vertex (Isolated_vertex *iv, Vertex* v)
  {
    iv->set_iterator(v);
    iso_verts.push_back(*iv);
  }

  /*! Erase an isloated vertex from the face. */
  void erase_isolated_vertex (Isolated_vertex *iv)
  {
    iso_verts.remove(*iv);
  }
};

/*! \class
 * Representation of an inner CCB.
 */
template <class Traits_, class CMap>
class Arr_cmap_inner_ccb
{
public:

  typedef Traits_ Traits;
  typedef typename CMap::Dart_handle Dart_handle;

public:

  typedef Arr_cmap_face<Traits, CMap>          Face;
  typedef Arr_cmap_halfedge<Traits, CMap>      Halfedge;
  typedef typename Face::Inner_ccb_iterator    Inner_ccb_iterator;

private:

  Face                *fa;   // The face the contains the CCB in its interior.
  Inner_ccb_iterator  hit;   // The inner CCB identifier.
  Dart_handle         dh;

public:

  /*default constructor*/
  Arr_cmap_hole() :
    fa(NULL)
  {}

  /*! Get a halfedge along the component (const version). */
  const Halfedge* halfedge () const
  {
    return (*hit);
  }

  /*! Get a halfedge along the component (non-const version). */
  Halfedge* halfedge ()
  {
    return (*hit);
  }

  /*! Set a representative halfedge for the component. */
  void set_halfedge (Halfedge *he)
  {
    *hit = he;
    dh = he->dh;
    return;
  }

  /*! Get the incident face (non-const version). */
  Face* face()
  {
    return fa;
  }

  /*! Get the incident face (const version). */
  const Face face() const
  {
    return fa;
  }

  /*! Get the iterator (non-const version). */
  Inner_ccb_iterator iterator()
  {
    return hit;
  }

  /*! Get the iterator (const version). */
  Inner_ccb_iterator iterator () const
  {
    return hit;
  }

  /*! Set the incident face. */
  void set_face(Face* f)
  {
    fa = f;
  }

  /*! Set the inner CCB iterator. */
  void set_iterator (Inner_ccb_iterator it)
  {
    hit = it;
    //dh = it->dh;
  }
};

/*! \class
 * Representation of an isolated vertex.
 */
template <class Traits_, class CMap>
class Arr_cmap_isolated_vertex
{
public:
    typedef Traits_ Traits;

public:
  typedef Arr_cmap_isolated_vertex<Traits, CMap>    Self;
  typedef Arr_cmap_face<Traits, CMap>               Face;
  typedef typename Face::Isolated_vertex_iterator   Isolated_vertex_iterator;

private:
  Face                      *f;      // The containing face.
  Isolated_vertex_iterator  ivit;    // The isolated vertex identifier.


public:
  /*! Default constructor. */
  Arr_cmap_isolated_vertex ():
    f(NULL)
  {}

  /*! Get the containing face (const version). */
  const Face* face () const
  {
    return (f);
  }

  /*! Get the containing face (non-const version). */
  Face* face ()
  {
    return (f);
  }

  /*! Set the incident face, the one that contains the isolated vertex. */
  void set_face (Face* fa)
  {
    f = fa;
  }

  /*! Get the isolated vertex iterator (const version). */
  Isolated_vertex_iterator iterator () const
  {
    //CGAL_assertion(iter_is_not_singular);
    return (ivit);
  }

  /*! Get the isolated vertex iterator (non-const version). */
  Isolated_vertex_iterator iterator ()
  {
    //CGAL_assertion(iter_is_not_singular);
    return (ivit);
  }

  /*! Set the isolated vertex iterator. */
  void set_iterator (Isolated_vertex_iterator iv)
  {
    ivit = iv;
    //iter_is_not_singular = true;
    //return;
  }
};

/*! \class
 * Representation of an outer CCB.
 */
template <class Traits_, class CMap>
class Arr_cmap_outer_ccb
{
public:

  typedef Traits_ Traits;

public:

  typedef Arr_cmap_outer_ccb<Traits, CMap>     Self;
  typedef Arr_cmap_halfedge<Traits,CMap>       Halfedge;
  typedef Arr_cmap_face<Traits, CMap>          Face;
  typedef typename Face::Outer_ccb_iterator    Outer_ccb_iterator;

private:

  Face               *f;    // The face the contains the CCB in its interior.
  Outer_ccb_iterator oit;   // The outer CCB identifier.
  Dart_handle        dh;

public:

  /*default constructor*/
  Arr_cmap_outer_ccb () {};

  /*! Get a halfedge along the component (non-const version). */
  Halfedge* halfedge()
  {
    return *oit;
  }

  /*! Get a halfedge along the component (const version). */
  const Halfedge* halfedge() const
  {
    return *oit;
  }

  /*! Set a representative halfedge for the component. */
  void set_halfedge (Halfedge *he)
  {
    *oit = he;
    dh = he->dh;
  }

  /*! Get the incident face (non-const version). */
  Face* face ()
  {
    return f;
  }

  /*! Get the incident face (const version). */
  const Face* face () const
  {
    return f;
  }

  /*! Set the incident face. */
  void set_face(Face* fa)
  {
    f = fa;
  }

  /*! Get the iterator (const version). */
  Outer_ccb_iterator iterator () const
  {
    return oit;
  }

  /*! Get the iterator (non-const version). */
  Outer_ccb_iterator iterator ()
  {
    //CGAL_assertion(iter_is_not_singular);
    return (oit);
  }

  /*! Set the outer CCB iterator. */
  void set_iterator (Outer_ccb_iterator it)
  {
    oit = it;
  }
};

template <class Traits_, class CMap>
class Arr_cmap_dcel
{
public:

  typedef Traits_                                    Traits;
  typedef Arr_cmap_dcel<Traits, CMap>                Self;
  typedef Arr_cmap_vertex<Traits, CMap>              Vertex;
  typedef Arr_cmap_halfedge<Traits, CMap>            Halfedge;
  typedef Arr_cmap_face<Traits, CMap>                Face;
  typedef Arr_cmap_inner_ccb<Traits, CMap>           Inner_ccb;
  typedef Arr_cmap_outer_ccb<Traits, CMap>           Outer_ccb;
  typedef Arr_cmap_isolated_vertex<Traits, CMap>     Isolated_vertex;
  typedef typename CMap::Dart_handle                 Dart_handle;

  typedef Inner_ccb                                  Hole;

public:

  CMap cm;

protected:

  typedef std::list<Vertex>              Vertex_list;
  typedef std::list<Halfedge>            Halfedge_list;
  typedef std::list<Face>                Face_list;
  typedef std::list<Hole>                Inner_ccb_list;
  typedef std::list<Outer_ccb>           Outer_ccb_list;
  typedef std::list<Isolated_Vertex>     Iso_vert_list;

public:

  typedef Halfedge_list::size_type Size;

protected:

  Vertex_list         vertices;     // The vertices container.
  Halfedge_list       halfedges;    // The halfedges container.
  Face_list           faces;        // The faces container.
  Outer_ccb_list      out_ccbs;     // The outer CCBs.
  Inner_ccb_list      in_ccbs;      // The inner CCBs.
  Iso_vert_list       iso_verts;    // The isolated vertices.

public:

  // Definitions of iterators.
  typedef typename Vertex_list::iterator               Vertex_iterator;
  typedef typename Halfedge_list::iterator             Halfedge_iterator;
  typedef typename Face_list::iterator                 Face_iterator;
  typedef CGAL::N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;

  // Definitions of const iterators.
  typedef typename Vertex_list::const_iterator        Vertex_const_iterator;
  typedef typename Halfedge_list::const_iterator      Halfedge_const_iterator;
  typedef typename Face_list::const_iterator          Face_const_iterator;
  typedef CGAL::N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

public:

  /*! Default constructor. */
  Arr_cmap_dcel()
  {}

  /*! Destructor. */
  ~Arr_cmap_dcel ()
  {
    delete_all();
  }

  /*! Get the number of DCEL vertices. */
  Size size_of_vertices () const
  {
    return (vertices.size());
  }

  /*! Get the number of DCEL halfedges (twice the number of edges). */
  Size size_of_halfedges () const
  {
    return (halfedges.size());
  }

  /*! Get the number of DCEL faces. */
  Size size_of_faces() const
  {
    return (faces.size());
  }

  /*! Get the number of outer CCBs. */
  Size size_of_outer_ccbs() const
  {
    return (out_ccbs.size());
  }

  /*! Get the number of inner CCBs. */
  Size size_of_inner_ccbs() const
  {
    return (in_ccbs.size());
  }

  /*! Get the number of isolated vertices. */
  Size size_of_isolated_vertices () const
  {
    return (iso_verts.size());
  }

  Vertex_iterator   vertices_begin()  { return vertices.begin(); }
  Vertex_iterator   vertices_end()    { return vertices.end(); }
  Halfedge_iterator halfedges_begin() { return halfedges.begin();}
  Halfedge_iterator halfedges_end()   { return halfedges.end(); }
  Face_iterator     faces_begin()     { return faces.begin(); }
  Face_iterator     faces_end()       { return faces.end(); }
  Edge_iterator     edges_begin()     { return halfedges.begin(); }
  Edge_iterator     edges_end()       { return halfedges.end(); }
    
  Vertex_const_iterator   vertices_begin() const { return vertices.begin(); }
  Vertex_const_iterator   vertices_end() const { return vertices.end(); }
  Halfedge_const_iterator halfedges_begin() const { return halfedges.begin(); }
  Halfedge_const_iterator halfedges_end() const { return halfedges.end(); }
  Face_const_iterator     faces_begin() const { return faces.begin(); }
  Face_const_iterator     faces_end() const { return faces.end(); }
  Edge_const_iterator     edges_begin() const { return halfedges.begin(); }
  Edge_const_iterator     edges_end() const { return halfedges.end(); }

  /*! Create a new vertex. */
  Vertex* new_vertex()
  {
    Vertex* v;
    //v->dcel = this;
    vertices.push_back(*v);
    return v;
  }

  /*! Create a new pair of opposite halfedges. */
  Halfedge* new_edge()
  {
    Halfedge* he1;
    Halfedge* he2;
    he1->dcel = this;
    he2->dcel = this;

    he1->dh = cm.create_dart();
    he2->dh = cm.create_dart();

    he1->set_opposite(he2);
    he2->set_opposite(he1);
    
    cm.template link_beta<2>(he1->dh, he2->dh);
    halfedges.push_back(*he1);
    halfedges.push_back(*he2);

    return he1;
  }

  /*! Create a new face. */
  Face* new_face()
  {
    Face *face;
    face->dcel = this;
    face->f = cm.create_face();
    faces.push_back(*face);
    return face;
  }

  /*! Create a new outer CCB. */
  Outer_ccb* new_outer_ccb ()
  {
    Outer_ccb *oc;
    out_ccbs.push_back(*oc);
    return oc;
  }

  /*! Create a new inner CCB. */
  Inner_ccb* new_inner_ccb ()
  {
    Inner_ccb *ic;
    in_ccbs.push_back(*ic);
    return ic;
  }

  /*! Create a new isolated vertex. */
  Isolated_vertex* new_isolated_vertex ()
  {
    Isolated_vertex *iv;
    iso_verts.push_back(*iv);
    return iv;
  }

  /*! Delete an existing vertex. */
  void delete_vertex (Vertex *v)
  {
    vertices.remove (*v);
    
  }

  /*! Delete an existing pair of opposite halfedges. */
  void delete_edge (Halfedge *h)
  {
    Halfedge *h_opp = h->opposite();
    delete_halfedge(*h);
    delete_halfedge(*h_opp);
    
  }

  /*! Delete an existing face. */
  void delete_face(Face *f)
  {
    faces.remove (*f);
    
  }

  /*! Delete an existing outer CCB. */
  void delete_outer_ccb (Outer_ccb *oc)
  {
    out_ccbs.remove (*oc);
    
  }

  /*! Delete an existing inner CCB. */
  void delete_inner_ccb (Inner_ccb *ic)
  {
    in_ccbs.remove (*ic);
    
  }

  void delete_halfedge (Halfedge *h)
  {
    halfedges.remove(*h);
    
  }

  /*! Delete an existing isolated vertex. */
  void delete_isolated_vertex (Isolated_vertex *iv)
  {
    iso_verts.remove(*iv);
    
  }

  /*! Delete all DCEL features. */
  void delete_all()
  {
    // Free all vertices.
    Vertex_iterator    vit = vertices.begin(), v_curr;

    while (vit != vertices.end())
    {
      v_curr = vit;
      ++vit;
      delete_vertex (&(*v_curr));
    }

    // Free all halfedges.
    Halfedge_iterator  hit = halfedges.begin(), h_curr;

    while (hit != halfedges.end())
    {
      h_curr = hit;
      ++hit;
      delete_halfedge (&(*h_curr));
    }
      
    // Free all faces.
    Face_iterator      fit = faces.begin(), f_curr;

    while (fit != faces.end())
    {
      f_curr = fit;
      ++fit;
      delete_face (&(*f_curr));
    }

    // Free all outer CCBs.
    typename Outer_ccb_list::iterator   ocit = out_ccbs.begin(), oc_curr;

    while (ocit != out_ccbs.end())
    {
      oc_curr = ocit;
      ++ocit;
      delete_outer_ccb (&(*oc_curr));
    }

    // Free all inner CCBs.
    typename Inner_ccb_list::iterator   icit = in_ccbs.begin(), ic_curr;

    while (icit != in_ccbs.end())
    {
      ic_curr = icit;
      ++icit;
      delete_inner_ccb (&(*ic_curr));
    }

    // Free all isolated vertices.
    typename Iso_vert_list::iterator   ivit = iso_verts.begin(), iv_curr;

    while (ivit != iso_verts.end())
    {
      iv_curr = ivit;
      ++ivit;
      delete_isolated_vertex (&(*iv_curr));
    }
  }

  /*!
   * Assign our DCEL the contents of another DCEL.
   */
  void assign (const Self& dcel)
  {
    
  }
};

}//namespace CGAL

#endif
