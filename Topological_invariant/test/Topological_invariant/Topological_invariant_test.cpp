#include <iostream>
#include <CGAL/Topological_surface.h>
#include <CGAL/exceptions.h>
#include <cassert>

template<class Map, class Dart_const_handle>
void printAlpha(const Map& map, Dart_const_handle dh, int i) {
    if(!map.is_free(dh, i)){
        Dart_const_handle d = map.alpha(dh, i);
        std::cout << "\t\talpha "<<i<<" : " << &*d << '\n';
    }
}

template<class Map>
void printMap(const Map& map) {
    typedef typename Map::Dart_const_handle Dart_const_handle;
    typedef typename Map::Dart_range Dart_range;
    typedef typename Dart_range::const_iterator Dart_const_iterator;

    for (Dart_const_iterator it = map.darts().begin(); it != map.darts().end(); ++it) {
        std::cout << "\tdart " << &*it << " "<< (map.dart_signature(it)?"+":"-") << std::endl;
        printAlpha(map, it, 0);
        printAlpha(map, it, 1);
        printAlpha(map, it, 2);
        
    }
}

void test_make_combinatorial_tetrahedron(){
    CGAL::Topological_surface<> surface;
    
    surface.make_combinatorial_tetrahedron();
    
    assert(surface.is_valid());
    std::cout<<"test_make_combinatorial_tetrahedron : OK"<<std::endl;
}

void test_make_combinatorial_hexahedron(){
    CGAL::Topological_surface<> surface;
    
    surface.make_combinatorial_hexahedron();
    
    assert(surface.is_valid());
    std::cout<<"test_make_combinatorial_hexahedron : OK"<<std::endl;
}

void test_make_combinatorial_polygon(){
    CGAL::Topological_surface<> surface;
    
    surface.make_combinatorial_polygon(3);
    
    assert(surface.is_valid());
    std::cout<<"test_make_combinatorial_polygon : OK"<<std::endl;
}

void test_vertex_next() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;


    Surface surface;
    
    surface.make_combinatorial_tetrahedron();
    
    Dart_handle d = surface.darts().begin();
    Halfedge_handle he = surface.halfedge(d);
    
    assert(surface.vertex_next(surface.vertex_next(surface.vertex_next(he)))==he);
    
    std::cout<<"test_vertex_next : OK"<<std::endl;
}

void test_link_alpha() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;

    Surface surface;

    Dart_handle a1 = surface.create_dart();
    Dart_handle a2 = surface.create_dart();
    Dart_handle b1 = surface.create_dart();
    Dart_handle b2 = surface.create_dart();

    surface.link_alpha<2>(a1, a2);
    surface.link_alpha<2>(b1, b2);

    surface.link_alpha<0>(a1, b1);
    surface.link_alpha<0>(a2, b2);

    assert(surface.is_valid());

    Halfedge_handle he_a1 = surface.halfedge(a1);
    Halfedge_handle he_a2 = surface.halfedge(a2);
    Halfedge_handle he_b1 = surface.halfedge(b1);
    Halfedge_handle he_b2 = surface.halfedge(b2);

    assert(he_a1 != Halfedge_handle());
    assert(he_a2 != Halfedge_handle());
    assert(he_b1 != Halfedge_handle());
    assert(he_b2 != Halfedge_handle());

    assert(he_a1 == he_a2);
    assert(he_b1 == he_b2);

    assert(surface.halfedge_opposite(he_a1) == he_b1);

    std::cout << "test_link_alpha : OK" << std::endl;
}

void testSurfaceNoLink0() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;

    Surface surface;

    /*std::cout<<"begin error"<<std::endl;
    Dart_handle a1 = surface.create_dart();

    try {
        assert(surface.halfedge_opposite(surface.halfedge(a1)) == Halfedge_handle());
    } catch (CGAL::Precondition_exception& e) {

    }
    std::cout<<"end error"<<std::endl;*/

    std::cout << "No link 0 : OK" << std::endl;
}

void testSurfaceNoLink1() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;

    Surface surface;

    Dart_handle a = surface.create_dart();
    Dart_handle b = surface.create_dart();
    Dart_handle c = surface.create_dart();

    //  c 0  b  1  a
    // ===|=====0===== 2

    surface.link_alpha<1>(a, b);
    surface.link_alpha<0>(b, c);

    Halfedge_handle he = surface.halfedge(a);

    assert(he != Halfedge_handle());

    std::cout << "No link 1 : OK" << std::endl;
}

void testSurfaceNoLink2() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;

    Surface surface;

    Dart_handle a = surface.create_dart();
    Dart_handle b = surface.create_dart();

    Halfedge_handle he_a = surface.halfedge(a);
    Halfedge_handle he_b = surface.halfedge(b);

    assert(he_a != Halfedge_handle());
    assert(he_b != Halfedge_handle());

    surface.link_alpha<2>(a, b);
    surface.unlink_alpha<2>(a);

    he_a = surface.halfedge(a);
    he_b = surface.halfedge(b);

    assert(he_a != Halfedge_handle());
    assert(he_b != Halfedge_handle());

    std::cout << "No Link 2 : OK" << std::endl;
}

void testHalfedgeFaceLoop() {
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;

    Surface surface;
    
    surface.make_combinatorial_polygon(3);
    
    
    assert(surface.is_valid());

    Halfedge_handle first = surface.halfedge(surface.darts().begin());
    
    if(surface.dart_signature(surface.dart(first))){
        first = surface.halfedge_opposite(first);
    }
    
    Halfedge_handle he = first;
    
    assert(!surface.dart_signature(surface.dart(he)));
    
    assert(surface.halfedge_signature(he));
    
    he = surface.halfedge_next(he);
    
    assert(!surface.dart_signature(surface.dart(he)));
    
    assert(surface.halfedge_signature(he));
    
    he = surface.halfedge_next(he);
    
    assert(!surface.dart_signature(surface.dart(he)));
    
    assert(surface.halfedge_signature(he));
    
    he = surface.halfedge_next(he);
    
    assert(first==he);

    std::cout << "halfedge face loop : OK" << std::endl;
}

void test_path(){
    typedef CGAL::Topological_surface<> Surface;
    typedef Surface::Dart_handle Dart_handle;
    typedef Surface::Halfedge_handle Halfedge_handle;
    typedef Surface::Path_handle Path_handle;

    Surface surface;
    
    surface.make_combinatorial_tetrahedron();
    
    Path_handle path = surface.create_path();
    
    Halfedge_handle he = surface.halfedge(surface.darts().begin());
    
    path->push_front(he);
    
    he = surface.vertex_next(he);
    
    path->push_front(he);

    std::cout << "test_path : OK" << std::endl;
    
}

int main(){
    test_make_combinatorial_tetrahedron();
    test_make_combinatorial_hexahedron();
    test_make_combinatorial_polygon();
    test_vertex_next();
    test_link_alpha();
    testSurfaceNoLink0();
    testSurfaceNoLink1();
    testSurfaceNoLink2();
    testHalfedgeFaceLoop();
    test_path();
    
    std::cout<<"Success"<<std::endl;
    return 0;
}
