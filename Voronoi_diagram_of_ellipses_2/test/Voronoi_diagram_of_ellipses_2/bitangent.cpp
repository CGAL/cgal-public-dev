// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>

#include <CGAL/Bitangent.h>

#include <CGAL/Timer.h>

typedef CGAL::Gmpfi INT;
typedef CGAL::Gmpq QQ;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;

Ellipse_2 e1, e2, e3;

using namespace std;

void test_bitangent()
{
    CGAL::VORELL::Bitangent<ET> bt12(e1, e2);
    cout << "external bitangent of {e1,e2}: " << ET::to_interval(bt12.external()) << endl;
    cout << "external bitangent of {e1,e2}: " << ET::to_interval(bt12.other_external()) << endl;
    cout << "internal bitangent of {e1,e2}: " << ET::to_interval(bt12.internal()) << endl;
    cout << "internal bitangent of {e1,e2}: " << ET::to_interval(bt12.other_internal()) << endl;
    
    cout << "ext/int bitangents of {e1,e2}: "  << CGAL::to_double(bt12.external()) << ", "
                                                << CGAL::to_double(bt12.internal()) << ", " 
                                                << CGAL::to_double(bt12.other_internal()) << ", " 
                                                << CGAL::to_double(bt12.other_external()) << endl;

    cout << "relative position of {e1,e2}: " << bt12.relative_position() << endl;
    cout << "degenerate relative position of {e1,e2}: " << bt12.is_degenerate_pair() << endl;
    
    CGAL::VORELL::Bitangent<ET> bt13(e1, e3);
    cout << "ext/int bitangents of {e1,e3}: "  << CGAL::to_double(bt13.external()) << ", "
                                                << CGAL::to_double(bt13.internal()) << ", " 
                                                << CGAL::to_double(bt13.other_internal()) << ", " 
                                                << CGAL::to_double(bt13.other_external()) << endl;
    cout << "relative position of {e1,e3}: " << bt13.relative_position() << endl;
    cout << "degenerate relative position of {e1,e3}: " << bt13.is_degenerate_pair() << endl;

    cout << "relative position of tanpoint internal(e1,e2) with e3: "
         << e3.boundary_relpos(e1, bt12.internal()) << endl;
    cout << "relative position of tanpoint internal(e1,e3) with e2: "
         << e2.boundary_relpos(e1, bt13.internal()) << endl;
}

int main(int argc, char **argv)
{
    std::cin >> e1 >> e2 >> e3;

//    CGAL::set_pretty_mode (cout);
    test_bitangent();

    return 0;
}

