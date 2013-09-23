
#ifndef files_Arr_cmap_dcel_h
#define files_Arr_cmap_dcel_h
/*! \file
 * The definition of the DCEL class based on Combinatorial_map_with_holes.
 */

#include <CGAL/Combinatorial_map_with_holes.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/Cell_attribute_with_point.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Dart.h>
#include <CGAL/Compact_container.h>
#include <list>
#include <iostream>

namespace CGAL {
/*! \class
 * The arrangement DCEL class represented by Combinatorial_map_with_holes.
 * The Traits parameters corresponds to a geometric traits class, which
 * defines the Point_2 and X_monotone_curve_2 types.
 */

/* Forward declarations. */
template <int d, class Refs> class Arr_cmap_halfedge;
template <class LCC,class Refs> class Arr_cmap_vertex;
template <unsigned int d_, class Traits_, class Items_, class Alloc_>
    class Arr_cmap_dcel;
template <unsigned int d_, class Refs> class Arr_cmap_outer_ccb;
template <unsigned int d_, class Refs> class Arr_cmap_inner_ccb;
template <unsigned int d_, class Refs> class Arr_cmap_isolated_vertex;

/* The min items for the Arr_cmap_dcel class. */
struct Arr_cmap_min_items
{
  /// Dart_wrapper defines the type of darts used, and enabled attributes.
  template < class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Arr_cmap_halfedge< 2, Refs > Dart;
    typedef CGAL::Cell_attribute<Refs, typename Refs::X_monotove_curve>
                                                                  Edge_attrib;
    typedef CGAL::cpp11::tuple<typename Refs::Vertex, Edge_attrib> Attributes;
  };
};

template <class LCC, class Refs>
class Arr_cmap_vertex : public Cell_attribute_with_point<LCC>
{
public:

  typedef Cell_attribute_with_point<LCC>           Base;
  typedef typename Base::Point                     Point;
  typedef typename Refs::Dart_handle               Dart_handle;
  typedef typename Refs::Dart_const_handle         Dart_const_handle;
  typedef typename Refs::Isolated_vertex           Isolated_vertex;

private:

  char pss[2];
  bool is_not_null;
  Isolated_vertex* iso;

public:

  /* default constructor */
  Arr_cmap_vertex () :
    is_not_null(false),
    iso(NULL)
  {
    pss[0] = pss[1] = static_cast<char> (CGAL::ARR_INTERIOR);
  }

  /*! Destructor. */
  virtual ~Arr_cmap_vertex() {}

  bool has_null_point () const
  {
    return (!is_not_null);
  }

  /*! Set the point (may be a NULL point). */
  void set_point (Point *p)
  {
    this->mpoint = (*p);
  }

  bool is_isolated () const
  {
    if (this->dart() == NULL)
      return true;
    else
      return false;
  }

  Dart_handle halfedge ()
  {
    CGAL_precondition (! is_isolated());
    return this->dart();
  }

  Dart_const_handle halfedge () const
  {
    CGAL_precondition (! is_isolated());
    return this->dart();
  }

  void set_halfedge (Dart_handle dh)
  {
    this->set_dart(dh);
  }

  /*! Get the boundary type in x. */
  Arr_parameter_space parameter_space_in_x () const
  {
    return (Arr_parameter_space (pss[0]));
  }

  /*! Get the boundary type in y. */
  Arr_parameter_space parameter_space_in_y () const
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
  virtual void assign (const Arr_cmap_vertex<LCC, Refs>& v)
  {
    is_not_null = v.is_not_null;
    this->mpoint = v.point();
    pss[0] = v.pss[0];
    pss[1] = v.pss[1];
  }

  /*! Get the isolated vertex information (const version). */
  const Isolated_vertex* isolated_vertex () const
  {
    CGAL_precondition (is_isolated());
    return iso;
  }

  /*! Get the isolated vertex information (non-const version). */
  Isolated_vertex* isolated_vertex ()
  {
    CGAL_precondition (is_isolated());
    return iso;
  }

  /*! Set the isolated vertex information. */
  void set_isolated_vertex (Isolated_vertex* iv)
  {
    iso = iv;
  }
};

template <int d, class Refs>
class Arr_cmap_halfedge: public Dart<d, Refs>
{
public:

  typedef Dart<d, Refs>                                    Base;
  typedef Arr_cmap_halfedge<d, Refs>                       Halfedge;
  typedef typename Refs::LCC                               LCC;
  typedef Arr_cmap_vertex<LCC, Refs>                       Vertex;
  typedef typename Refs::template Attribute_handle<0>::type
                                                           Vertex_handle;
  typedef typename Refs::template Attribute_const_handle<0>::type
                                                           Vertex_const_handle;
  typedef typename Refs::Dart_handle                       Dart_handle;
  typedef typename Refs::Dart_const_handle                 Dart_const_handle;
  typedef typename Refs::X_monotove_curve                  X_monotone_curve;

  typedef typename Refs::Outer_ccb                         Outer_ccb;
  typedef typename Refs::Inner_ccb                         Inner_ccb;

private:

  bool left_to_right;    // Determine the direction of the halfedge
  bool on_inner_ccb;     // Determine whether it is on inner ccb
  bool is_not_null;
  Outer_ccb  *outer;
  Inner_ccb *inner;

public:

  bool has_null_curve () const
  {
    return (!is_not_null);
  }

  /*! Get the x-monotone curve (const version). */
  const X_monotone_curve& curve() const
  {
    CGAL_precondition (!is_not_null);
    return (this->template attribute<1>()->info());
  }

  /*! Get the x-monotone curve (non-const version). */
  X_monotone_curve& curve ()
  {
    CGAL_precondition (!is_not_null);
    return (this->template attribute<1>()->info());
  }

  /*! Set the x-monotone curve. */
  void set_curve (X_monotone_curve* c)
  {
    this->template attribute<1>()->info() = (*c);
    is_not_null = true;
  }

  /*! Get the previous halfedge along the chain (const version). */
  Dart_const_handle prev () const
  {
    CGAL_assertion(!is_free(0));
    return this->beta(0);
  }

  /*! Get the previous halfedge along the chain (non-const version). */
  Dart_handle prev ()
  {
    CGAL_assertion(!is_free(0));
    return this->beta(0);
  }

  /*! Get the next halfedge along the chain (const version). */
  Dart_const_handle next () const
  {
    CGAL_assertion(!is_free(1));
    return this->beta(1);
  }

  /*! Get the next halfedge along the chain (non-const version). */
  Dart_handle next ()
  {
    CGAL_assertion(!is_free(1));
    return this->beta(1);
  }

  /*! Sets the opposite halfedge. */
  void set_opposite (Dart_handle dh)
  {
    this->template basic_link_beta<2>(dh);
  }

  /*! Set the previous halfedge along the chain. */
  void set_prev (Dart_handle dh)
  {
    this->template basic_link_beta<0>(dh);
  }

  /*! Set the next halfedge along the chain. */
  void set_next (Dart_handle dh)
  {
    this->template basic_link_beta<1>(dh);
  }

  Vertex_handle vertex ()
  {
    return this->template attribute<0>();
  }

  Vertex_const_handle vertex () const
  {
    return this->template attribute<0>();
  }

  void set_vertex(Vertex* v)
  {
    this->template attribute<0>()->info() = (*v);
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
  void set_direction(Arr_halfedge_direction dir)
  {
    if (dir == ARR_LEFT_TO_RIGHT)
    {
      left_to_right = true;
      this->beta(2)->left_to_right = false;
    }
    else
    {
      left_to_right = false;
      this->beta(2)->left_to_right = true;
    }
  }

  /*! Assign from another halfedge. */
  virtual void assign (const Arr_cmap_halfedge<d, Refs>& he)
  {
    is_not_null = he.is_not_null;
    if (is_not_null)
    {
      this->template attribute<1>()->info() = he.template attribute<1>->info();
    }
  }

  /*! Check whether the halfedge lies on the boundary of an inner CCB. */
  bool is_on_inner_ccb () const
  {
    return on_inner_ccb;
  }

  /*!
   * Get an incident outer CCB (const version).
   * \pre The edge does not lie on an inner CCB.
   */
  const Outer_ccb* outer_ccb () const
  {
    CGAL_precondition (! is_on_inner_ccb());
    return outer;
  }

  /*!
   * Get an incident outer CCB (non-const version).
   * \pre The edge does not lie on an inner CCB.
   */
  Outer_ccb* outer_ccb ()
  {
    CGAL_precondition (! is_on_inner_ccb());
    return outer;
  }

  /*! Set the incident outer CCB. */
  void set_outer_ccb (Outer_ccb *oc)
  {
    outer = oc;
    on_inner_ccb = false;
  }

  /*!
   * Get an incident inner CCB (const version).
   * \pre The edge lies on an inner CCB.
   */
  const Inner_ccb* inner_ccb () const
  {
    CGAL_precondition (is_on_inner_ccb());
    return inner;
  }

  /*!
   * Get an incident inner CCB (non-const version).
   * \pre The edge lies on an inner CCB.
   */
  Inner_ccb* inner_ccb ()
  {
    CGAL_precondition (is_on_inner_ccb());
    return inner;
  }

  /*! Set the incident inner CCB. */
  void set_inner_ccb (Inner_ccb *ic)
  {
    inner = ic;
    on_inner_ccb = true;
  }
};

template <unsigned int d_, class Refs>
class Arr_cmap_outer_ccb
{
public:

  typedef typename Refs::Face                        Face;
  typedef typename Face::Outer_ccb_iterator          Outer_ccb_iterator;
  typedef typename Refs::Face_iterator               Face_iterator;
  typedef typename Refs::Dart_handle                 Dart_handle;
  typedef typename Refs::Dart_const_handle           Dart_const_handle;

private:

  Face_iterator           fit;
  Outer_ccb_iterator      iter;
  bool                    iter_is_not_singular;

public:

  /* default constructor */
  Arr_cmap_outer_ccb () :
    iter_is_not_singular(false)
  {}

  /*! Copy constructor. */
  Arr_cmap_outer_ccb (const Arr_cmap_outer_ccb& other )
    : iter_is_not_singular(other.iter_is_not_singular)
  {
    fit = other.fit;
    if(other.iter_is_not_singular)
    {
      iter = other.iter;
    }
  }

  /*! Get a halfedge along the component (const version). */
  Dart_handle halfedge () const
  {
    return (iter);
  }

  /*! Get a halfedge along the component (non-const version). */
  Dart_handle halfedge ()
  {
    return (iter);
  }

  /*! Set a representative halfedge for the component. */
  void set_halfedge (Dart_handle he)
  {
    iter = he;
    return;
  }

  /*! Get the incident face (const version). */
  Face_iterator face () const
  {
    return (fit);
  }

  /*! Get the incident face (non-const version). */
  Face_iterator face ()
  {
    return (fit);
  }

  /*! Set the incident face. */
  void set_face (Face_iterator f)
  {
    fit = f;
    return;
  }

  /*! Get the iterator (const version). */
  Outer_ccb_iterator iterator () const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Get the iterator (non-const version). */
  Outer_ccb_iterator iterator ()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Set the outer CCB iterator. */
  void set_iterator (Outer_ccb_iterator it)
  {
    iter = it;
    iter_is_not_singular = true;
    return;
  }
};

template <unsigned int d_, class Refs>
class Arr_cmap_inner_ccb
{
public:

  typedef typename Refs::Face                        Face;
  typedef typename Face::Inner_ccb_iterator          Inner_ccb_iterator;
  typedef typename Refs::Face_iterator               Face_iterator;
  typedef typename Refs::Face_const_iterator         Face_const_iterator;
  typedef typename Refs::Dart_handle                 Dart_handle;
  typedef typename Refs::Dart_const_handle           Dart_const_handle;
    
private:
    
  Face_iterator           fit;
  Inner_ccb_iterator      iter;
  bool                    iter_is_not_singular;

public:

  /* default constructor */
  Arr_cmap_inner_ccb () :
    iter_is_not_singular(false)
  {}

  /*! Copy constructor. */
  Arr_cmap_inner_ccb (const Arr_cmap_inner_ccb& other )
    : iter_is_not_singular(other.iter_is_not_singular)
  {
    fit = other.fit;
    if(other.iter_is_not_singular)
    {
      iter = other.iter;
    }
  }

  /*! Get a halfedge along the component (const version). */
  Dart_handle halfedge () const
  {
    return (iter);
  }

  /*! Get a halfedge along the component (non-const version). */
  Dart_handle halfedge ()
  {
    return (iter);
  }

  /*! Set a representative halfedge for the component. */
  void set_halfedge (Dart_handle he)
  {
    iter = he;
    return;
  }

  /*! Get the incident face (const version). */
  Face_iterator face () const
  {
    return (fit);
  }

  /*! Get the incident face (non-const version). */
  Face_iterator face ()
  {
    return (fit);
  }

  /*! Set the incident face. */
  void set_face (Face_iterator f)
  {
    fit = f;
    return;
  }

  /*! Get the iterator (const version). */
  Inner_ccb_iterator iterator () const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Get the iterator (non-const version). */
  Inner_ccb_iterator iterator ()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Set the outer CCB iterator. */
  void set_iterator (Inner_ccb_iterator it)
  {
    iter = it;
    iter_is_not_singular = true;
    return;
  }
};

template <unsigned int d_, class Refs>
class Arr_cmap_isolated_vertex
{
public:

  typedef typename Refs::Face                        Face;
  typedef typename Face::Isolated_vertex_iterator
                                                     Isolated_vertex_iterator;
  typedef typename Refs::Face_iterator               Face_iterator;

private:

  Face_iterator              fit;
  Isolated_vertex_iterator   iv_it;
  bool                       iter_is_not_singular;

public:

  /* default constructor */
  Arr_cmap_isolated_vertex () :
    iter_is_not_singular(false)
  {}

  /*! Copy constructor. */
  Arr_cmap_isolated_vertex (const Arr_cmap_isolated_vertex& other )
    : iter_is_not_singular(other.iter_is_not_singular)
  {
    fit = other.fit;
    if(other.iter_is_not_singular)
    {
      iv_it = other.iv_it;
    }
  }

  /*! Get the containing face (const version). */
  Face_iterator face () const
  {
    return (fit);
  }

  /*! Get the containing face (non-const version). */
  Face_iterator face ()
  {
    return (fit);
  }

  /*! Set the incident face, the one that contains the isolated vertex. */
  void set_face (Face_iterator f)
  {
    fit = f;
    return;
  }

  /*! Get the isolated vertex iterator (const version). */
  Isolated_vertex_iterator iterator () const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iv_it);
  }

  /*! Get the isolated vertex iterator (non-const version). */
  Isolated_vertex_iterator iterator ()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iv_it);
  }

  /*! Set the isolated vertex iterator. */
  void set_iterator (Isolated_vertex_iterator iv)
  {
    iv_it = iv;
    iter_is_not_singular = true;
    return;
  }
};

template <unsigned int d_, class Traits_,
          class Items_=Arr_cmap_min_items,
          class Alloc_=CGAL_ALLOCATOR(int)>
class Arr_cmap_dcel: public Combinatorial_map_with_holes<d_, Items_, Alloc_>
{
public:

  typedef Traits_                                            Traits;
  typedef Items_                                             Items;
  typedef Alloc_                                             Alloc;

  typedef typename Traits::X_monotone_curve_2                X_monotone_curve;
  typedef typename Traits::Kernel                            Kernel;

  typedef Combinatorial_map_with_holes<d_, Items, Alloc>     Base;
  typedef Arr_cmap_dcel<d_, Traits, Items, Alloc>            Refs;

  typedef Linear_cell_complex_traits<2, Kernel>              LCC_Traits;
  typedef Linear_cell_complex<2, 2, LCC_Traits>              LCC;
  typedef typename LCC::Point                                Point;

  typedef Arr_cmap_vertex<LCC, Refs>                         Vertex;
  typedef Arr_cmap_halfedge<d_, Refs>                        Halfedge;
  typedef Arr_cmap_outer_ccb<d_, Refs>                       Outer_ccb;
  typedef Arr_cmap_inner_ccb<d_, Refs>                       Inner_ccb;
  typedef Arr_cmap_isolated_vertex<d_, Refs>                 Isolated_vertex;

  typedef Inner_ccb                                          Hole;

  typedef typename Base::Dart_handle                         Dart_handle;
  typedef typename Base::Dart_const_handle                   Dart_const_handle;

 // typedef typename Base::Dart_container                      Dart_container;
  typedef typename Items::template Dart_wrapper<Refs>        Dart_wrapper;
  typedef typename Dart_wrapper::Dart                        Dart;
  typedef typename Alloc::template rebind<Dart>::other       Dart_allocator;
  typedef Compact_container<Dart,Dart_allocator>             Dart_container;

  typedef typename Base::Face_iterator                       Face_iterator;
  typedef typename Base::Face_const_iterator                 Face_const_iterator;
  typedef CGAL::N_step_adaptor_derived<Dart_handle, 2>       Edge_iterator;
  typedef CGAL::N_step_adaptor_derived<Dart_const_handle, 2> Edge_const_iterator;

  typedef typename Dart_container::size_type                 Size;
  typedef typename Dart_container::size_type                 size_type;
  typedef typename Dart_container::difference_type           difference_type;
  typedef typename Dart_container::difference_type           Difference;

  typedef typename Refs::template Attribute_range<0>::type   Vertex_list;
  typedef typename Vertex_list::iterator                     Vertex_iterator;
  typedef typename Vertex_list::const_iterator          Vertex_const_iterator;
  typedef typename Refs::template Attribute_handle<0>::type  Vertex_handle;

  typedef Compact_container<Outer_ccb>                       Outer_ccb_list;
  typedef typename Outer_ccb_list::iterator                  Out_ccb_iterator;
  typedef typename Outer_ccb_list::const_iterator            Out_ccb_const_iterator;

  typedef Compact_container<Inner_ccb>                       Inner_ccb_list;
  typedef typename Inner_ccb_list::iterator                  Inn_ccb_iterator;
  typedef typename Inner_ccb_list::const_iterator            Inn_ccb_const_iterator;

  typedef Compact_container<Isolated_vertex>                 Iso_vert_list;
  typedef typename Iso_vert_list::iterator                   Iso_vert_iterator;
  typedef typename Iso_vert_list::const_iterator             Iso_vert_const_iterator;

protected:

  //Vertex_list                 vertices;      // The vertices container.
  Outer_ccb_list              out_ccbs;      // The outer CCBs.
  Inner_ccb_list              in_ccbs;       // The inner CCBs.
  Iso_vert_list               iso_verts;     // The isolated vertices.

public:

  class Face: public Base::Face
  {
  public:

    typedef Dart_container                       Outer_ccbs_container;
    typedef typename Outer_ccbs_container::iterator
                                                 Outer_ccb_iterator;
    typedef typename Outer_ccbs_container::const_iterator
                                                 Outer_ccb_const_iterator;

    typedef Dart_container                       Inner_ccbs_container;
    typedef typename Inner_ccbs_container::iterator
                                                 Inner_ccb_iterator;
    typedef typename Inner_ccbs_container::const_iterator
                                                 Inner_ccb_const_iterator;

    typedef Vertex_list                          Isolated_vertices_container;
    typedef typename Isolated_vertices_container::iterator
                                                 Isolated_vertex_iterator;
    typedef typename Isolated_vertices_container::cons_iterator
                                                 Isolated_vertex_const_iterator;

  protected:

    enum
    {
      IS_UNBOUNDED = 1,
      IS_FICTITIOUS = 2
    };

    int                            flags;      // Face flags.
    Outer_ccbs_container           outer_ccbs; // The outer CCBs of the faces.
    Inner_ccbs_container           inner_ccbs; // The inner CCBs of the face.
    Isolated_vertices_container    iso_verts;  // The isolated vertices inside
                                               // the face.

  public:

    /*default constructor*/
    Face() :
      flags(0)
    {}

    /*! Check if the face is unbounded. */
    bool is_unbounded () const
    {
      return ((flags & IS_UNBOUNDED) != 0);
    }

    /*! Set the face as bounded or unbounded. */
    void set_unbounded (bool unbounded)
    {
      flags = (unbounded) ? (flags | IS_UNBOUNDED) : (flags & ~IS_UNBOUNDED);
    }

    /*! Check if the face is fictitious. */
    bool is_fictitious () const
    {
      return ((flags & IS_FICTITIOUS) != 0);
    }

    /*! Set the face as fictitious or valid. */
    void set_fictitious (bool fictitious)
    {
      flags = (fictitious) ? (flags | IS_FICTITIOUS) : (flags & ~IS_FICTITIOUS);
    }

    /*! Assign from another face. */
    virtual void assign (const Face& f)
    {
      flags = f.flags;
    }

  public:

    typedef Arr_cmap_vertex<LCC, Refs>              Vertex;
    typedef Arr_cmap_halfedge<d_, Refs>             Halfedge;
    typedef Arr_cmap_outer_ccb<d_, Refs>            Outer_ccb;
    typedef Arr_cmap_inner_ccb<d_, Refs>            Inner_ccb;
    typedef Arr_cmap_isolated_vertex<d_, Refs>      Isolated_vertex;

    typedef Inner_ccb                               Hole;

    /*! Get the number of outer CCBs the face has. */
    size_t number_of_outer_ccbs () const
    {
      return (this->outer_ccbs.size());
    }

    /*! Get an iterator for the first outer CCB of the face. */
    Outer_ccb_iterator outer_ccbs_begin()
    {
      return (this->outer_ccbs.begin());
    }

    /*! Get a past-the-end iterator for the outer CCBs inside the face. */
    Outer_ccb_iterator outer_ccbs_end()
    {
      return (this->outer_ccbs.end());
    }

    /*! Get an const iterator for the first outer CCB inside the face. */
    Outer_ccb_const_iterator outer_ccbs_begin() const
    {
      return (this->outer_ccbs.begin());
    }

    /*! Get a const past-the-end iterator for the outer CCBs inside the face. */
    Outer_ccb_const_iterator outer_ccbs_end() const
    {
      return (this->outer_ccbs.end());
    }
      
    /*! Add an outer CCB to the face. */
    void add_outer_ccb (Outer_ccb *oc, Dart_handle h)
    {
      Dart temp = *h;
      oc->set_iterator (this->outer_ccbs.insert (temp));
      return;
    }

    /*! Erase an outer CCB of the face. */
    void erase_outer_ccb (Outer_ccb *oc)
    {
      this->outer_ccbs.erase (oc->iterator());
    }

    typedef Inner_ccb_iterator                        Hole_iterator;
    typedef Inner_ccb_const_iterator                  Hole_const_iterator;

    /*! Get the number of inner CCBs the face has. */
    size_t number_of_inner_ccbs () const
    {
      return (this->inner_ccbs.size());
    }

    /*! Get an iterator for the first inner CCB of the face. */
    Inner_ccb_iterator inner_ccbs_begin()
    {
      return (this->inner_ccbs.begin());
    }

    /*! Get a past-the-end iterator for the inner CCBs inside the face. */
    Inner_ccb_iterator inner_ccbs_end()
    {
      return (this->inner_ccbs.end());
    }

    /*! Get an const iterator for the first inner CCB inside the face. */
    Inner_ccb_const_iterator inner_ccbs_begin() const
    {
      return (this->inner_ccbs.begin());
    }

    /*! Get a const past-the-end iterator for the inner CCBs inside the face. */
    Inner_ccb_const_iterator inner_ccbs_end() const
    {
      return (this->inner_ccbs.end());
    }

    /*! Add an inner CCB to the face. */
    void add_inner_ccb (Inner_ccb *ic, Dart_handle h)
    {
      Dart temp = *h;
      ic->set_iterator (this->inner_ccbs.insert (temp));
      return;
    }

    /*! Erase an inner CCB of the face. */
    void erase_inner_ccb (Inner_ccb *ic)
    {
      this->inner_ccbs.erase (ic->iterator());
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
      return (this->iso_verts.size());
    }

    /*! Get an iterator for the first isloated vertex inside the face. */
    Isolated_vertex_iterator isolated_vertices_begin()
    {
      return (this->iso_verts.begin());
    }

    /*! Get a past-the-end iterator for the isloated vertices inside the face. */
    Isolated_vertex_iterator isolated_vertices_end()
    {
      return (this->iso_verts.end());
    }

    /*! Get an const iterator for the first isloated vertex inside the face. */
    Isolated_vertex_const_iterator isolated_vertices_begin() const
    {
      return (this->iso_verts.begin());
    }

    /*! Get a const past-the-end iterator for the isloated vertices inside the
     * face. */
    Isolated_vertex_const_iterator isolated_vertices_end() const
    {
      return (this->iso_verts.end());
    }

    /*! Add an isloated vertex inside the face. */
    void add_isolated_vertex (Isolated_vertex *iv, Vertex_iterator v)
    {
        
      iv->set_iterator (this->iso_verts.insert (*v));
      return;
    }

    /*! Erase an isloated vertex from the face. */
    void erase_isolated_vertex (Isolated_vertex *iv)
    {
      this->iso_verts.erase (iv->iterator());
      return;
    }
  };

public:

  //Vertex allocator
  typedef typename Alloc::template rebind<Vertex>        Vertex_alloc_rebind;
  typedef typename Vertex_alloc_rebind::other            Vertex_allocator;

  // Face allocator.
  typedef typename Alloc::template rebind<Face>          Face_alloc_rebind;
  typedef typename Face_alloc_rebind::other              Face_allocator;

  // Outer CCB allocator.
  typedef typename Alloc::template rebind<Outer_ccb>     Out_ccb_alloc_rebind;
  typedef typename Out_ccb_alloc_rebind::other           Outer_ccb_allocator;

  // Inner CCB allocator.
  typedef typename Alloc::template rebind<Inner_ccb>     In_ccb_alloc_rebind;
  typedef typename In_ccb_alloc_rebind::other            Inner_ccb_allocator;

  // Isolated vertex allocator.
  typedef typename Alloc::template rebind<Isolated_vertex>
                                                         Iso_vert_alloc_rebind;
  typedef typename Iso_vert_alloc_rebind::other          Iso_vert_allocator;

protected:

  Outer_ccb_allocator out_ccb_alloc;        // An allocator for outer CCBs.
  Inner_ccb_allocator in_ccb_alloc;         // An allocator for inner CCBs.
  Iso_vert_allocator  iso_vert_alloc;       // Allocator for isolated vertices.

private:

  // Copy constructor - not supported.
  Arr_cmap_dcel (const Refs&);

  // Assignment operator - not supported.
  Refs& operator = (const Refs&);

public:

  /*! Default constructor. */
  Arr_cmap_dcel ()
  {}

  /*! Destructor. */
  ~Arr_cmap_dcel ()
  {
    delete_all();
  }

  /*! Get the number of DCEL vertices. */
  Size size_of_vertices () const
  {
    return this->template number_of_attributes<0>();
  }

  /*! Get the number of DCEL halfedges (twice the number of edges). */
  Size size_of_halfedges () const
  {
    return this->number_of_darts();
  }

  /*! Get the number of DCEL faces. */
  Size size_of_faces() const
  {
    return this->facets.size();
  }

  /*! Get the number of outer CCBs. */
  Size size_of_outer_ccbs() const
  {
    return out_ccbs.size();
  }

  /*! Get the number of inner CCBs. */
  Size size_of_inner_ccbs() const
  {
    return in_ccbs.size();
  }

  /*! Get the number of isolated vertices. */
  Size size_of_isolated_vertices () const
  {
    return iso_verts.size();
  }

public:

  Vertex_iterator   vertices_begin()
  { return this->template attribute<0>().begin(); }
  Vertex_iterator   vertices_begin()
  { return this->template attribute<0>().end() }
  Face_iterator     faces_begin()         { return this->facets.begin(); }
  Face_iterator     faces_end()           { return this->facets.end(); }
  Dart_handle       halfedges_begin()     { return this->first_dart(); }
  Dart_handle       halfedges_end()       { return this->mdart.end(); }
  Edge_iterator     edges_begin()         { return this->first_dart(); }
  Edge_iterator     edges_end()           { return this->mdart.end(); }


  Vertex_const_iterator vertices_begin() const
  { return this->template attribute<0>().begin(); }
  Vertex_const_iterator vertices_end() const
  { return this->template attribute<0>().end() }
  Dart_const_handle     halfedges_begin() const { return this->first_dart(); }
  Dart_const_handle     halfedges_end() const { return this->mdart.end(); }
  Face_const_iterator   faces_begin() const { return this->facets.begin(); }
  Face_const_iterator   faces_end() const { return this->facets.end(); }
  Edge_const_iterator   edges_begin() const { return this->first_dart(); }
  Edge_const_iterator   edges_end() const { return this->mdart.end(); }

public:

  /*! Create a new vertex. */
  Vertex_handle new_vertex()
  {
    Vertex_handle vh = this->template create_attribute<0>();
    return vh;
  }

  /*! Create a new pair of opposite halfedges. */
  Dart_handle new_edge()
  {
    Dart_handle h1 = this->create_dart();
    Dart_handle h2 = this->create_dart();
/*
    this->template set_attribute<0>(h1, this->template create_attribute<0>());
    this->template set_attribute<0>(h2, this->template create_attribute<0>());
*/
    this->template set_attribute<1>(h1, this->template create_attribute<1>());
    this->template set_attribute<1>(h2, this->template create_attribute<1>());

    h1->set_opposite(h2);
    h2->set_opposite(h1);
    return h1;
  }

  /*! Create a new face. */
  Face_iterator new_face()
  {
    Face_iterator fit = this->create_face();
    return fit;
  }

  Outer_ccb* new_outer_ccb ()
  {
    Outer_ccb  *oc = out_ccb_alloc.allocate (1);
    out_ccb_alloc.construct (oc, Outer_ccb());
    oc = &*out_ccbs.insert (*oc);
    return (oc);
  }

  /*! Create a new inner CCB. */
  Inner_ccb* new_inner_ccb ()
  {
    Inner_ccb  *ic = in_ccb_alloc.allocate (1);
    in_ccb_alloc.construct (ic, Inner_ccb());
    ic = &*in_ccbs.insert (*ic);
    return (ic);
  }

  /*! Create a new isolated vertex. */
  Isolated_vertex* new_isolated_vertex ()
  {
    Isolated_vertex  *iv = iso_vert_alloc.allocate (1);
    
    iso_vert_alloc.construct (iv, Isolated_vertex());
    iv = &*iso_verts.insert (*iv);
    return (iv);
  }

  void delete_vertex (Vertex_handle v)
  {
    this->template erase_attribute<0>(v);
  }

  void delete_edge (Dart_handle dh)
  {
    Dart_handle dh2 = dh->beta(2);
    this->erase_dart(dh);
    this->erase_dart(dh2);
  }

  void delete_face(Face_iterator f)
  {
    this->delete_face(f);
  }

  void delete_outer_ccb (Outer_ccb_iterator oc)
  {
    out_ccbs.erase (oc);
  }

  void delete_inner_ccb (Inner_ccb_iterator ic)
  {
    in_ccbs.erase (ic);
  }

  void delete_isolated_vertex (Isolated_vertex_iterator iv)
  {
    iso_verts.erase (iv);
  }

  void delete_all()
  {
    out_ccbs.clear();
    in_ccbs.clear();
    iso_verts.clear();
    this->facets.clear();
    this->clear();
  }
};

}//namespace CGAL

#endif
