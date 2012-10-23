// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef GENERATOR_DEFS_H
#define GENERATOR_DEFS_H

template < class NT >
struct Generator {  // represents a generator matrix (GJG^T or GB^T)
    Generator() {
    }

    Generator(unsigned n) :
        a(n), b(n), c(n), d(n) {
    }
        
    //! copies \c n elements from respective iterators to [a,b,c,d]
    template < class InputIterator >
    Generator(unsigned n, InputIterator ai, InputIterator bi, 
            InputIterator ci, InputIterator di) :
            a(n), b(n), c(n), d(n) {

        _copy_data(ai, bi, ci, di, 0, n);
    }

    template < class InputIterator >
    Generator(unsigned sz, unsigned n, InputIterator ai, InputIterator bi,
            InputIterator ci, InputIterator di) :
            a(sz, NT(0)), b(sz, NT(0)), c(sz, NT(0)), d(sz, NT(0)) {

        _copy_data(ai, bi, ci, di, 0, n);
    }

    //! copies the range of \c [fe;be) from \c g to \c this
    //! it is assumed that the destination has enough space for copy
    void assign_range(unsigned fe, unsigned be, const Generator& g) {

        _copy_data(g.a.begin(), g.b.begin(), g.c.begin(), g.d.begin(), fe, be);
    }

    // shifts all vectors to the left by the amount given by n
    void shift_left(unsigned n) {
        unsigned s = a.size() - n;
//         if((int)s <= 0)

        _copy_data(a.begin(), b.begin(), c.begin(), d.begin(),
                s, a.size());
        a.resize(s); b.resize(s); c.resize(s); d.resize(s);
    }

protected:

    // copies the range [first; beyond] from (ai,bi,ci,di) to *this
    template < class InputIterator >
    void _copy_data(InputIterator ai, InputIterator bi, InputIterator ci,
                InputIterator di, unsigned first, unsigned beyond) {

        std::copy(ai + first, ai + beyond, a.begin());
        std::copy(bi + first, bi + beyond, b.begin());
        std::copy(ci + first, ci + beyond, c.begin());
        std::copy(di + first, di + beyond, d.begin());
    }

public:
    std::vector< NT > a, b, c, d;

    void print() const {
        printf("\n========== G(a, b): ============\n");
        print_vector(a);
        print_vector(b);
        printf("========== B(c, d): ============\n");
        print_vector(c);
        print_vector(d);
    }
};

template < class NT >
struct GB_update { // updates for non-symmetric generator pair

    GB_update() {
    }
    GB_update(unsigned n) : a(n), c(n) {
    }
    std::vector< NT > a, c;
};

template < class NT >
struct GJG_update :
        public std::vector< NT > {

    typedef std::vector< NT > Base;

    GJG_update() : Base() {
    }
    GJG_update(unsigned n) : Base(n) {
    }
    GJG_update(unsigned n, const NT& x) : Base(n, x) {
    }
};

#endif // GENERATOR_DEFS_H
