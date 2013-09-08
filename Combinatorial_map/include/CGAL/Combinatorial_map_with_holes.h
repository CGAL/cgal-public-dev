
#ifndef files_Combinatorial_map_with_holes_h
#define files_Combinatorial_map_with_holes_h

#include <iostream>
#include <list>
#include <CGAL/Combinatorial_map.h>
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

  typedef Combinatorial_map_with_holes<d_, Items_,Alloc_>  Self;
  typedef CGAL::Combinatorial_map_base<d_, Self, Items_, Alloc_> Base;
  typedef typename Base::Dart_handle Dart_handle;
  typedef typename Base::Dart_const_handle Dart_const_handle;
  typedef typename Base::Alloc Alloc;

public:
  class Face
  {
  public:

    std::list<Dart_handle>    dart_list;

  public:

    /* default constructor */
    Face () {}

    /* Add the associated outer_boundary to the face */
    void add_boundary(Dart_handle dh) {
      if (dart_list.size() == 0) dart_list.push_back(dh);
      else dart_list.front() = dh;
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
      dart_list.front() = Base::null_dart_handle;
    }

    /* Erase the hole of the face */
    void erase_hole(Dart_handle dh)
    {
      dart_list.remove(dh);
    }
  };

public:

  std::list<Face>    facets;
      
public:

  /* Define an instance of Face, and return this Face object */
  Face create_face() {
    Face fa;
    facets.push_back(fa);
    return fa;
  }
};
  
}//namespace

#endif
