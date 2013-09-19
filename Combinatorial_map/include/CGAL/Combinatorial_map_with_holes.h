
#ifndef files_Combinatorial_map_with_holes_h
#define files_Combinatorial_map_with_holes_h

#include <iostream>
#include <list>
#include <vector>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Compact_container.h>
namespace CGAL {

template < unsigned int d_, 
           class Items_=CGAL::Combinatorial_map_min_items<d_>,
           class Alloc_=CGAL_ALLOCATOR(int) >
class Combinatorial_map_with_holes :
    public CGAL::Combinatorial_map_base<d_,
                            Combinatorial_map_with_holes<d_,Items_,Alloc_>,
                            Items_, Alloc_>
{
public:

  typedef Combinatorial_map_with_holes<d_, Items_,Alloc_>   Self;
  typedef CGAL::Combinatorial_map_base<d_, Self, Items_, Alloc_>
                                                            Base;
  typedef typename Base::Dart_handle                        Dart_handle;
  typedef typename Base::Dart_const_handle                  Dart_const_handle;
  typedef typename Base::Alloc                              Alloc;

public:

  class Face
  {

  public:

    std::vector<Dart_handle>    dart_list;

  public:

    /* default constructor */
    Face () {}

    /* Add the associated outer_boundary to the face */
    void add_boundary(Dart_handle dh) {
      if (dart_list.size() == 0) dart_list.push_back(dh);
      else dart_list[0] = dh;
    }

    /* Add a hole into the face */
    void add_hole(Dart_handle dh)
    {
      if (dart_list.size() == 0) {
          dart_list.push_back(Base::null_dart_handle);
        dart_list.push_back(dh);
      }
      else dart_list.push_back(dh);
    }

    /* Erase the outer_boundary of the face */
    void erase_boundary()
    {
      if (dart_list.size() == 1)
          dart_list.clear();
      else dart_list[0] = Base::null_dart_handle;
    }

    /* Erase the hole of the face */
    void erase_hole(Dart_handle dh)
    {
      for (int i = 0; i<dart_list.size(); ++i)
      {
        if (dart_list[i] == dh)
        {
          dart_list.erase(dh);
          break;
        }
      }
    }

  public:

    void * for_compact_container() const
    { return vp; }
    void * & for_compact_container()
    { return vp; }

  private:
    void * vp;

  };

public:
    
    typedef Compact_container<Face>                      Face_container;
    typedef typename Face_container::iterator            Face_handle;
    typedef typename Face_container::const_iterator      Face_const_handle;
    typedef typename Face_container::size_type           size_type;

public:

  Face_container       facets;

public:

  /* Define an instance of Face, and return this Face object */
  Face_handle create_face() {
    Face fa;
    return facets.template emplace<Face>(fa);
  }

  void delete_face(Face_handle fit)
  {
    facets.erase(fit);
  }

  size_type number_of_faces() const
  {
    return facets.size();
  }
};
  
}//namespace

#endif
