// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_POINT_GENERATORS_3_H
#define CGAL_POINT_GENERATORS_3_H 1
#include <CGAL/generators.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/internal/weighted_random_element.h>
#include <vector>
#include <algorithm>

#include <iostream>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

namespace CGAL {

template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef P result_type;
    typedef Random_points_in_sphere_3<P,Creator> This;
    Random_points_in_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open sphere with radius r, i.e. |`*g'| < r .
        // Three random numbers are needed from `rnd' for each point
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    r *= std::pow( this->_rnd.get_double() , 1.0/3.0 );  
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_on_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef P result_type;
    typedef Random_points_on_sphere_3<P,Creator> This;
    Random_points_on_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the sphere with radius r, i.e. |`*g'| == r . A
        // two random numbers are needed from `rnd' for each point.
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_on_sphere_3<P,Creator>::
generate_point() {
  // A strip between z and z+dz has an area independant of z
    typedef typename Creator::argument_type T;
    double alpha = this->_rnd.get_double() * 2.0 * CGAL_PI;
    double z     = 2 * this->_rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    Creator creator;
    this->d_item = creator( T(this->d_range * r * std::cos(alpha)),
                            T(this->d_range * r * std::sin(alpha)),
                            T(this->d_range * z));
}


template < class P, class Creator = 
                   Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_cube_3 : public Random_generator_base<P>{
    void generate_point();
public:
    typedef P result_type;
    typedef Random_points_in_cube_3<P,Creator> This;
    Random_points_in_cube_3( double a = 1, Random& rnd = default_random)
    : Random_generator_base<P>( a, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_cube_3<P,Creator>::
generate_point() {
    typedef typename Creator::argument_type T;
    Creator creator;
    this->d_item =
	     creator( T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)),
                      T(this->d_range * ( 2 * this->_rnd.get_double() - 1.0)));
}


template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, 
                         OutputIterator o, Creator creator)
{
    if  (n == 0)
        return o;

    int m = int(std::ceil(
                  std::sqrt(std::sqrt(static_cast<double>(n)))));

    while (m*m*m < int(n)) m++;

    double base = -a;  // Left and bottom boundary.
    double step = 2*a/(m-1);
    int j = 0;
    int k = 0;
    double px = base;
    double py = base;
    double pz = base;
    *o++ = creator( px, py, pz);
    for (std::size_t i = 1; i < n; i++) {
        j++;
        if ( j == m) {
           k++;
           if ( k == m) {
              py = base;
              px = base;
              pz = pz + step;
              k = 0;
           }
           else {
              px = base;
              py = py + step;
           }
           j = 0;
        } else {
           px = px + step;
        }
        *o++ = creator( px, py, pz);
    }
    return o;
}

template <class OutputIterator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, OutputIterator o)
{
    typedef std::iterator_traits<OutputIterator> ITraits;
    typedef typename ITraits::value_type         P;
    return points_on_square_grid_3(a, n, o, 
                 Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P>());
}

template < class P, class Creator = 
Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_triangle_3 : public Random_generator_base<P> {
	P _p,_q,_r;
	void generate_point();
public:
	typedef P result_type;
	typedef Random_points_in_triangle_3<P> This;
	typedef typename Kernel_traits<P>::Kernel::Triangle_3 Triangle_3;
	Random_points_in_triangle_3() {}
	Random_points_in_triangle_3( const This& x,Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd ),_p(x._p),_q(x._q),_r(x._r) {
		generate_point();
	}
	Random_points_in_triangle_3( const P& p, const P& q, const P& r, Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd ),_p(p),_q(q),_r(r) {
		generate_point();
	}
	Random_points_in_triangle_3( const Triangle_3& triangle,Random& rnd = default_random)
	: Random_generator_base<P>( 1,
			rnd),_p(triangle[0]),_q(triangle[1]),_r(triangle[2]) {
		generate_point();
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};


template<class P, class Creator >
void Random_points_in_triangle_3<P, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	Creator creator;
	double a1 = this->_rnd.get_double(0,1);
	double a2 = this->_rnd.get_double(0,1);
	if(a1>a2) std::swap(a1,a2);
	double b1 = a1;
	double b2 = a2-a1;
	double b3 = 1.0-a2;
	T ret[3];
	for(int i = 0; i < 3; ++i) {
	    ret[i] = T(to_double(_p[i])*b1+to_double(_q[i])*b2+to_double(_r[i])*b3);
	}
	this->d_item = creator(ret[0],ret[1],ret[2]);
}

template < class P, class Creator = 
Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_tetrahedron_3 : public Random_generator_base<P> {
	P _p,_q,_r,_s;
	void generate_point();
public:
	typedef P result_type;
	typedef Random_points_in_tetrahedron_3<P> This;
	typedef typename Kernel_traits<P>::Kernel::Tetrahedron_3 Tetrahedron_3;
	Random_points_in_tetrahedron_3() {}
	Random_points_in_tetrahedron_3( const This& x,Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd ),_p(x._p),_q(x._q),_r(x._r),_s(x._s) {
		generate_point();
	}
	Random_points_in_tetrahedron_3( const P& p, const P& q, const P& r, const P& s,Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd ),_p(p),_q(q),_r(r),_s(s) {
		generate_point();
	}
	Random_points_in_tetrahedron_3( const Tetrahedron_3& tetrahedron,Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd),_p(tetrahedron[0]),_q(tetrahedron[1]),_r(tetrahedron[2]),_s(tetrahedron[3]) {
		generate_point();
	}
	This operator=(This x) {
		_p = x._p;
		_q = x._q;
		_r = x._r;
		_s = x._s;
		return *this;
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};

template<class P, class Creator >
void Random_points_in_tetrahedron_3<P, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	Creator creator;
	double a[3];
	for(int i = 0; i < 3; ++i) {
		a[i]=this->_rnd.get_double(0,1);
	}
	std::sort(a,a+3);
	double b[4];
	b[0]=a[0];
	b[1]=a[1]-a[0];
	b[2]=a[2]-a[1];
	b[3]=1.0-a[2];
	T ret[3];
	for(int i = 0; i < 3; ++i) {
	    ret[i] = T(to_double(_p[i])*b[0]+to_double(_q[i])*b[1]+to_double(_r[i])*b[2]+to_double(_s[i])*b[3]);
	}
	this->d_item = creator(ret[0],ret[1],ret[2]);
}

/*
// @param:
// 	n = number of points that will be generated
template<typename OutputIterator, typename Element_RandomAccessIterator,
	typename VolumeElementFunctor, typename PointGeneratorFunctor>
OutputIterator
Random_points_in_mesh_3( Element_RandomAccessIterator el_begin,
		Element_RandomAccessIterator el_end,)
{
	int N = el_end - el_begin;
	std::vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> > container;
	container.reserve(N);
	Element_RandomAccessIterator it = el_begin;
	int i = 0;
	for (; it != el_end; it++) {
		VolumeElementFunctor volElem(*it);
		double weight = volElem();
		double presum = (i == 0 ? weight : weight +
				container[i-1].getPresum());
		PointGeneratorFunctor randGen(*it);
		CGAL::internal::Weighted_random_element<PointGeneratorFunctor> aux(randGen, presum);
		container[i] = aux;
		i++;
	}

	CGAL::Random rand;
	for (int i = 0; i < n; i++) {
		double tmp_presum = rand.get_double(0,
				container[N-1].getPresum());
		CGAL::internal::Weighted_random_element<PointGeneratorFunctor>
			tmp(tmp_presum);
		typename std::vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> >::iterator SampleIterator = upper_bound(container.begin(), container.end(), tmp);
		int SampleIndex = SampleIterator - container.begin();
		Element_RandomAccessIterator SampleElement = el_begin + SampleIndex;
		// *o++ = container[SampleIndex].getRand()();
		CGAL::cpp11::copy_n(container[SampleIndex].getRand,1, std::back_inserter(o));
	}
	return o;
}
*/

// @template
// 	P = type of points
// 	C3t3 = type of mesh
template < class Tri, class P, class C3t3, class VolumeElementFunctor, class
PointGeneratorFunctor, class Creator = 
Creator_uniform_3<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_mesh_3 : public Random_generator_base<P> {
	C3t3 _c3t3;
	void generate_point();
public:
	typedef Random_points_in_mesh_3<Tri, P, C3t3, VolumeElementFunctor,
		PointGeneratorFunctor> This;
	Random_points_in_mesh_3() {}
	Random_points_in_mesh_3( const C3t3& c3t3, Random& rnd = default_random)
	: Random_generator_base<P>( 1, rnd ),_c3t3(c3t3) {
		generate_point();
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};


template<class Tri, class P, class C3t3, class VolumeElementFunctor, class
PointGeneratorFunctor, class Creator >
void Random_points_in_mesh_3<Tri, P,C3t3, VolumeElementFunctor,
     PointGeneratorFunctor, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	typedef typename Kernel_traits<P>::Kernel::Tetrahedron_3 Tetrahedron_3;
	Creator creator;

	std::cout << "Actual number of cells in _c3t3 in complex: " <<
		_c3t3.number_of_cells_in_complex() << "\n";
	
	Tri tr = _c3t3.triangulation();
	int Nr_cells = _c3t3.number_of_cells_in_complex();
	Tetrahedron_3 *tetra;
	tetra = new Tetrahedron_3[Nr_cells];
	int i = 0;
	typename Tri::Finite_cells_iterator iter = tr.finite_cells_begin();
	for ( ; iter != tr.finite_cells_end(); ++iter) {
		if (_c3t3.is_in_complex(iter)) {
			tetra[i] = tr.tetrahedron(iter);
			i++;
		} else {
			Nr_cells--;
		}
	}


	std::vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> > container;
	container.reserve(Nr_cells);
	for (int i = 0; i < Nr_cells; i++) {
		VolumeElementFunctor volElem(tetra[i]);
		double weight = volElem();
		double presum = (i == 0 ? weight : weight +
				container[i-1].getPresum());
		PointGeneratorFunctor randGen(tetra[i][0], tetra[i][1],
				tetra[i][2], tetra[i][3]);
		CGAL::internal::Weighted_random_element<PointGeneratorFunctor> aux(randGen, presum);
		container[i] = aux;
	}

	CGAL::Random rand;
	double tmp_presum = rand.get_double(0,
			container[Nr_cells-1].getPresum());
	CGAL::internal::Weighted_random_element<PointGeneratorFunctor>
		tmp(tmp_presum);
	typename std::vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> >::iterator SampleIterator = upper_bound(container.begin(), container.end(), tmp);
	int SampleIndex = SampleIterator - container.begin();
	Tetrahedron_3 SampleElement = tetra[SampleIndex];

	std::vector<P> ret;
	ret.reserve(1);

	CGAL::cpp11::copy_n(container[SampleIndex].getRand(),1,
			std::back_inserter(ret));

	this->d_item = creator(ret[0].x(),ret[1].y(),ret[2].z());
	delete[] tetra;
}

} //namespace CGAL

#endif // CGAL_POINT_GENERATORS_3_H //
// EOF //
