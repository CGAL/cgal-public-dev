// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
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
    
#ifndef _MATRIX_ALGEBRA_C_
#define _MATRIX_ALGEBRA_C_

#include <include/modular_arithm.h>

namespace CGAL {

namespace internal {

//! evaluates a polynomial f at \c x given by coefficients from an array \c U :
//! f[0][degx] f[1][degx] f[2][degx] .. f[degy][degx] (with \c m_stride )
//! f[0][degx-1] f[1][degx-1] f[2][degx-1] .. f[degy][degx-1]
//! ...
template < class OutputIterator >
OutputIterator poly_eval2_pure(const unsigned *U, zmod x, unsigned degx,
        unsigned degy, unsigned m_stride, OutputIterator oi) {

    unsigned i, j;
    const unsigned *pU, *ppU;
    
    for(i = 0, pU = U; i <= degy; i++, pU++) {
        zmod y(pU[0]);

        for(j = 1, ppU = pU + m_stride; j <= degx; j++, ppU += m_stride) {
            y = y * x + zmod(ppU[0]);
        }
        *oi++ = y;
    }
    return oi;
}

template < class NT >
void print_inv_vector(const std::vector< NT >& a,
        const NT& denom, unsigned i_start, unsigned i_beyond) {        
 
    unsigned n = i_beyond - i_start;
    std::vector< NT > aa(n);

    NT inv(mod_inverse(denom.x, zmod::modulus));
    
    for(unsigned i = i_start; (int)i < (int)i_beyond; i++) {
        aa[i - i_start] = a[i] * inv;
    }   
    print_vector(aa, true);
}

#if 0
//! returns pseudo-remainder
template < class NT >
void euclid_pseudo_division(const std::vector< NT >& A_,
        const std::vector< NT >& B_) {

    std::vector< NT > A = A_, B = B_;
    unsigned m = A.size()-1, n = B.size()-1, i, j, k = 0;
    unsigned quo_deg = m - n;

    std::vector< NT > R = A, Q(quo_deg + 1, NT(0));
    unsigned e = m - n + 1, r = m;
    
    NT d = B[n]; // leading coeff
    while(r >= n) { // r - degree of the divisor

        NT lcr = R[r];
        unsigned newr = -1u;
        // go over all B coeffs (r >= n); compute: R = R * d - S * B
        for(i = n, j = r; (int)j >= 0; i--, j--) {

            NT subs(0);
            if((int)i >= 0)
                subs = lcr * B[i];
                            
            R[j] = R[j] * d - subs;
            if(newr == -1u && R[j] != NT(0))
                newr = j;
        }
        if(newr == -1u) // remainder vanished
            newr = 0;

//         printf("\n%d: newr: %d: ", k, newr);
//         print_vector(R, false);

        for(i = 0; (int)i <= quo_deg; i++)
            Q[i] = Q[i] * d;
        Q[r - n] = lcr; // adds 'S' to Q
        r = newr, k++, e--;
        if(r == 0) break;
    }

    NT q = d.pow(e);
    for(i = 0; i <= quo_deg; i++) {
        Q[i] = Q[i] * q;
    }
    for(i = 0; i <= r; i++) {
        R[i] = R[i] * q;
    }
}
#endif

template < class NT >
unsigned pgcd_check(std::vector< NT >& A,
         std::vector< NT >& B, unsigned *out = 0,
            unsigned o_stride = 1) {

    std::vector< NT > *pR = &A, *pB = &B;
    unsigned m = A.size()-1, n = B.size()-1, r = m, i, j, k = 0,
            nv = (n + 4) & ~3; // 'nv' = (n + 1) 4-aligned controls the offset

    if(m < n) {
        std::swap(pR, pB);
        std::swap(r, n);
    }

    unsigned div_counter = 0, ss;
    while(1) {
        
    ss = 0;
    NT d = (*pB)[n], lcr; // leading coeff
    while(r >= n) { // r - degree of the divisor

        unsigned r2 = -1;
        lcr = (*pR)[r];
        // go over all B coeffs (r >= n); compute: R = R * d - S * B
        for(i = n, j = r; (int)j >= 0; i--, j--) {

            NT subs(0);
            if((int)i >= 0)
                subs = lcr * (*pB)[i];

            (*pR)[j] = (*pR)[j] * d - subs;
            if(r2 == -1u && (*pR)[j].x != 0)
                r2 = j;
        }

        ss++; //! NOTE we do not check for vanishing lcr here !!
        r = (r2 == -1u ? 0 : r2);

// printf("\nr: %d; n: %d; #iters: %d\n", r, n, ss);
// // // 
//         if(div_counter == 0 && ss == 1) {
//            printf("\nr: %d; n: %d; #iters: %d; next lcr: %d\n", r, n, ss,
//                 (*pR)[r].x);
// //            fillout(d, out, o_stride, r+1);
// //             unsigned ofs = (r >= nv ? 0 : (nv - r)*o_stride);
//             unsigned ofs = 0;
//             writeout((*pR), out + ofs, o_stride, r-n, r+1);
//             return 1;   
//         }

//         printf("\nr: %d; n: %d; #iters: %dx\n", r, n, ss);
        if(r == 0)
            break;
    }
        
//     printf("\nr: %d; n: %d; #iters: %d\n", r, n, ss);

    div_counter++;
    if(div_counter == 2000 || (
        r == 0 && (*pR)[0] == NT(0))) { // pB is a gcd
//         printf("###################### GCD found: %d; counter: %d; "
//             "next lcr/lcf: %d %d\n", 
//                 n, div_counter, (*pR)[r].x, (*pB)[n].x);
        writeout((*pB), out /*+ (nv - n)*o_stride*/, o_stride, 0, n+1);
//         writeout((*pR), out + 1*o_stride, o_stride, 0, r+1);
        break;
    }

//     if(div_counter >= 1) {
//        printf("\nr: %d; n: %d; #iters: %d\n", r, n, ss);
// //       fillout((*pR)[r], out, o_stride, r+1);
// //         unsigned ofs = (r >= nv ? 0 : (nv - r)*o_stride);
//         unsigned ofs = 0;        
//         writeout((*pR), out + ofs * o_stride, o_stride, 0, r+1);
// //          writeout((*pB), out + ofs * o_stride, o_stride, 0, n+1);
//          return 1;   
//     }
    
        std::swap(pR, pB); // swap the pointers and degrees
        std::swap(r, n);
    } // while(1)

//     writeout((*pB), out, o_stride, 0, n+1);

//     printf("\ndiv_counter: %d; total iters: %d\n", div_counter, ss);
    return (n + 1);
}

//! returns \c true if polynomials are coprime
//! NOTE: input vectors are modified by the algorithm
template < class NT >
bool euclid_pseudo_gcd(std::vector< NT >& A, std::vector< NT >& B) {

    std::vector< NT > *pR = &A, *pB = &B;
    unsigned m = A.size()-1, n = B.size()-1, r = m, i, j, k = 0;

    if(m < n) {
        std::swap(pR, pB);
        std::swap(r, n);
    }

    //! r - current degree of divisor
    //! n - current degree of dividend
    while(1) {

    NT d = (*pB)[n]; // leading coeff
    while(r >= n) { // r - degree of the divisor

        NT lcr = (*pR)[r];
        // in fact coefficient R[r] needs not be computed (because it vanishes)
        unsigned newr = -1u;
        // go over all B coeffs (r >= n); compute: R = R * d - S * B
//         for(i = n, j = r; (int)i >= 0; i--, j--) {
        for(i = n, j = r; (int)j >= 0; i--, j--) {

            NT subs(0);
            if((int)i >= 0)
                subs = lcr * (*pB)[i];
                            
            (*pR)[j] = (*pR)[j] * d - subs;
            if(newr == -1u && (*pR)[j] != NT(0))
                newr = j;
        }
        if(newr == -1u) // remainder vanished
            newr = 0;
 
        r = newr; //! r is a degree of remainder
        if(r == 0)
            break;
    }

//      printf("\n step: %d ", k++);
//        printf("\n current prem: ");
//         print_inv_vector(*pR, NT(1), 0, r+1);
        if(r == 0 && (*pR)[0] == NT(0)) // pB is a gcd
            break;

        std::swap(pR, pB); // swap the pointers and degrees
        std::swap(r, n);
    }
/*
    printf("\n============ gcd deg: %d\n", n);
    print_inv_vector(*pB, (*pB)[n], 0, n+1);*/
    return (n == 0);
}

//! computes monic resultant, then multiplies it by v1lc^m * v2lc^n
//! \c v1/2lc - leading coefficients of the polynomials
//! \c failed : indicates if the algorithm computes zero denominator
template < class NT >
NT resultant_modular_sylvester(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed,
            const limb *testit = NULL) {

    failed = false;
    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;

    zmod v1lc = v1[n];
    zmod inv1(mod_inverse(v1lc.x, zmod::modulus)), inv2(1);
    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));

    // [a, b] = G; [c, d] = H
    a[0] = NT(1);
    b[m] = NT(0)-NT(1);

    for(i = 0; i < r; i++) {
        c[i] = (i <= n ? v1[n - i] * inv1 : NT(0));
        d[i] = (i <= m ? v2[m - i] * inv2 : NT(0));
        if(i >= m)
            d[i] = d[i] - v1[r - i] * inv1;
    }

    NT la(1), lc(1), det(1), denom(1), v12lc_exp(1);
    j = 0;

#if 1
    for(j = 0; j < m; j++) {  // suppose m is even (at first)
        NT a0 = a[j], b0 = NT(0), c0 = c[j], d0 = d[j];

        v12lc_exp = v12lc_exp * v1lc;
        for(i = j + 1; i < r; i++) {
            d[i] = (d[i] - c[i]*d0);
        }

        a[m] = d0; // b[m] == -1
        for(i = r - 2; (int)i >= (int)j; i--) {
            a[i + 1] = a[i]; // shiftdown
            c[i + 1] = c[i];
        }
        b[j] = 0, d[j] = 0;
        a[j] = 0; c[j] = 0;
    }
#endif
//       printf("\n================== beginning iter: G:\n");
//       print_vector(a); print_vector(b);
// 
//       printf("\n B:\n");
//       print_vector(c); print_vector(d);

    for(; j < r; j++) {

    NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

    for(i = j; i < r; i++) {

        NT ai;
        ai = la*(a[i] * c0 - b[i]*d0);
        b[i] = lc*(-a[i]*b0 + b[i]*a0);
        a[i] = ai;

        NT ci = lc*(c[i]*a0 - d[i]*b0);
        d[i] = la*(-c[i]*d0 + d[i]*c0);
        c[i] = ci;
    }

    lc = la * lc * lc; la = a[j];
    det = det * c[j], denom = denom * lc;

    if(denom.x == 0) {
        printf("ill-conditioned matrix at j = %d\n", j);
        failed = true;
        return NT(0);
    }

    for(i = r - 2; (int)i >= (int)j; i--) {
        a[i + 1] = a[i]; // shiftdown
        c[i + 1] = c[i];
    }

    a[j] = 0; c[j] = 0;
    }

    det = det * v12lc_exp;

//     printf("modulus: %x det = %x\n", zmod::modulus, det.x);

#if CUMP_RUN_RESULTANTS_ONLY
    denom = zmod(1);
#endif
//!BUG BUG: weird results unless mod_inverse is declared as 'inline'
     NT inv(mod_inverse(denom.x, zmod::modulus));
     return NT(det * inv);
}

//! reconstructs a polynomial f(x) of degree \c xs.size() - 1 s.t.
//! \c f(xs[i]) = \c ys[i]
//! \c xs.size() == \c ys.size()
//! \c n - number of interpolation points: xs.size() >= n
//! \c stride - output data stride
template < class NT >
void vandermonde_interpolate(NT *xs, NT *ys, NT *res,
         unsigned n, unsigned stride = 1) {

    std::vector< zmod > a(2*n, zmod(0)), b(2*n, zmod(0));

    unsigned i, j;
    for(i = 0; i < n; i++) {
        a[i] = zmod(1); b[i] = zmod(ys[i]);
    }

    a[n] = zmod(1);
    zmod denom(1);

    for(j = 0; j < n; j++) { // first m steps are positive

    zmod a0 = a[j], b0 = b[j];
    denom = denom * a0;

    // 2*n length generator G, and n+1 length generator B
    for(i = j + 1; i < j + n; i++) {
        b[i] = -a[i] * b0 + b[i] * a0;
    }

    b[j + n] = -b0;

    zmod x0 = zmod(xs[j]);
    for(i = j + 1; i < n; i++) { // n times - however total is a lot...
    // xs[j] is x[1]
        a[i] = (zmod(xs[i]) - x0) * a[i];
    }
    a[j] = 0;

    // j = 0..n-1
    // a[n + j + 2] till a[2n - 1] are always zeros !!
    
    for(i = n + j; i >= n + 1; i--) { // n-1 times
        a[i] = a[i - 1] - x0 * a[i]; // to protect from rewriting
    }

    if(j < n - 1)
        a[n + j + 1] = NT(1); // a[n + j + 1] is constantly 1 !!
    // here i = n + j
    a[n] = -x0 * a[n];
    }

    // compute -denom^-1
    //! for the time being inv is not computed
#if CUMP_RUN_RES_AND_INTERPOLATE
    zmod inv(1);
#else
    zmod inv = -denom.pow(zmod::modulus - 2);
#endif

    //! should be ofs = j + 1 if you test intermediate results 
    unsigned ofs = n;  //! ofs = n at the end
    unsigned *pres = res;

    for(i = ofs; i < n + ofs; i++, pres += stride) {
        pres[0] = (b[i] * inv).x;
    }
}

//! \c res can be equal to \c xs or \c ys for in-place mod
inline void vandermonde_interpolate_poly(unsigned *xs, unsigned *ys,
        std::vector< zmod >& res, unsigned n) {

    std::vector< unsigned > rr(n);
    vandermonde_interpolate(xs, ys, rr.data(), n, 1);
    for(unsigned i = 0; i < n; i++)
        res[i] = zmod(rr[i]);
}

template < class NT >
NT sylvester_blockschur_gcd(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed,
            const limb *testit = NULL) {

    failed = false;
    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;

    NT v1lc = v1[n];
///!!!!     NT inv1(mod_inverse(v1lc.x, zmod::modulus)), inv2(1);
    NT inv1(1), inv2(1);
    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));
    std::vector< NT > ax(r, NT(0)), bx(r, NT(0)), cx(r, NT(0)), dx(r, NT(0));

    // [a, b] = G; [c, d] = B
    a[0] = NT(1);
    b[m] = NT(1);

    for(i = 0; i < r; i++) {
        c[i] = (i <= n ? v1[n - i] * inv1 : NT(0));
        d[i] = (i <= m ? v2[m - i] * inv2 : NT(0));
        if(i >= m)
            d[i] = d[i] - v1[r - i] * inv1;
    }

    print_vector(a);
    print_vector(b);
    print_vector(c);
    print_vector(d);

// NOTE: we do not need the block Schur algorithm during first m steps
// because the first m minors of Sylvester matrix are non-vanishing
// (determined by the leading coefficient of either polynomial)

    NT res(1), denom(1), ld(1); // ld - current denominator
    // NOTE: first make sure that r is even !!
    for(j = 0; j < r; j += 2) {
        // here: g: [x1 y1]   b: [z1 w1]
        //          [x2 y2]      [z2 w2]
        NT x1 = a[j], y1 = b[j], x2 = a[j + 1], y2 = b[j + 1];
        NT z1 = c[j], w1 = d[j], z2 = c[j + 1], w2 = d[j + 1];
        
        // compute G.b^+ and g.B^+ note that the leading 2x2 block
        // is same for both matrix products
        for(i = j; i < r; i++) {
            ax[i] = a[i] * z1 + b[i] * w1;
            bx[i] = a[i] * z2 + b[i] * w2;

            cx[i] = c[i] * x1 + d[i] * y1;
            dx[i] = c[i] * x2 + d[i] * y2;
        }

        NT Lx1 = ax[j], Ly1 = bx[j], Lx2 = ax[j + 1], Ly2 = bx[j + 1];
        NT Uz1 = cx[j], Uw1 = dx[j], Uz2 = cx[j + 1], Uw2 = dx[j + 1];

        // diagonal block D is given by: [Lx1    Ly1   ]
        //                               [Lx2 Ly2 + Lx1]
        res = res * (Lx1 * (Ly2 + Lx1) - Lx2 * Ly1);
        denom = denom * ld * ld * ld * ld;
 
        if(Lx1 != Uz1 || Ly1 != Uz2 || Lx2 != Uw1 || Ly2 != Uw2) {
            printf("fatal j = %d\n", j);
            printf("Lx1: %x; Ly1: %x; Lx2: %x; Ly2: %x\n", Lx1.x, Ly1.x,
                Lx2.x, Ly2.x);
            printf("Uz1: %x; Uw1: %x; Uz2: %x; Uw2: %x\n", Uz1.x, Uw1.x,
                Uz2.x, Uw2.x);
            throw 1;
        }

//     std::cout << "================== G.b^+:\n"; 
//     print_vector(ax);
//     print_vector(bx);
//     std::cout << "================== g.B^+:\n";
//     print_vector(cx);
//     print_vector(dx);

        // compute L2(ax,bx) and U2(cx,dx)^+, s.t., L2 - Z.L2.z^+ = G.b^+ and
        // U2 - z.U2.Z^+ = g.B^+ and multiply them by (In - Z) from the left
        for(i = r - 1; (int)i >= j + 1; i--) {
            ax[i] = ax[i] - ax[i - 1];
            cx[i] = cx[i] - cx[i - 1];
        }
        for(i = r - 1; (int)i >= j + 1; i--) {
            bx[i] = bx[i] - bx[i - 1] + ax[i - 1];
            dx[i] = dx[i] - dx[i - 1] + cx[i - 1];
        }

        Lx1 = ax[j], Ly1 = bx[j], Lx2 = ax[j + 1], Ly2 = bx[j + 1];
        Uz1 = cx[j], Uw1 = dx[j], Uz2 = cx[j + 1], Uw2 = dx[j + 1];
 
        // compute the inverse of leading 2x2 block of either matrix
        NT det = Lx1*Ly2 - Lx2*Ly1;
        // these two determinats are the same:
//         NT det1 = Lx1*Ly2 - Lx2*Ly1, det2 = Uz1*Uw2 - Uz2*Uw1; 
//         printf("det1: %x det2: %x\n", det1.x, det2.x);

        // compute  D_L^(-1).g where D_L = SubMatrix(L2, 2x2)
        NT gx1 = Ly2 * x1 - Ly1 * x2, gy1 = Ly2 * y1 - Ly1 * y2,
           gx2 = Lx1 * x2 - Lx2 * x1, gy2 = Lx1 * y2 - Lx2 * y1;
        // compute  D_U^(-1).b where D_U = SubMatrix(U2, 2x2)
        NT bz1 = Uw2 * z1 - Uw1 * z2, bw1 = Uw2 * w1 - Uw1 * w2,
           bz2 = Uz1 * z2 - Uz2 * z1, bw2 = Uz1 * w2 - Uz2 * w1;

        // subtract from original G and B 
        for(i = j; i < r; i++) {
            
            NT ai = (ax[i] * gx1 + bx[i] * gx2);
            bx[i] = (ax[i] * gy1 + bx[i] * gy2);
            ax[i] = ai;
    
            NT ci = (cx[i] * bz1 + dx[i] * bz2);
            dx[i] = (cx[i] * bw1 + dx[i] * bw2);
            cx[i] = ci;
            // now we just need to subtract these from original (G, B)
        }
        
        for(i = j; i < r; i++) {
            a[i] = a[i] * det - ax[i];
            b[i] = b[i] * det - bx[i];
            c[i] = c[i] * det - cx[i];
            d[i] = d[i] * det - dx[i];
        }
        ld = det * ld; // because of multiplication D^(-1).g

        NT invdet(mod_inverse(ld.x, zmod::modulus));
        for(i = j; i < r; i++) {
            ax[i] = a[i] * invdet;
            bx[i] = b[i] * invdet;
            cx[i] = c[i] * invdet;
            dx[i] = d[i] * invdet;
        }
//         std::cout << "\n";
//         for(unsigned k = j; k < r; k++) {
//             for(i = j; i < r; i++) {
//                 NT cc = ax[k] * cx[i] + bx[k] * dx[i];
//                 std::cout << cc << "  ";
//             }
//             std::cout << "|\n";
//         }
//      std::cout << " update G:\n"; 
//     print_vector(ax);
//     print_vector(bx);
//     std::cout << " update B:\n";
//     print_vector(cx);
//     print_vector(dx);

    }

    NT inv(mod_inverse(denom.x, zmod::modulus));
    res = res * inv;
    std::cout << "## resultant: " << res << "\n";

    return NT(0);
}

//! newer resultant version using simplified generator expressions
template < class NT >
NT resultant_modular_sylvester_new(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed,
            const limb *testit = NULL) {

    // polynomial degrees
    unsigned p = v1.size()-1, q = v2.size()-1, r = p + q, i, j;
    if(p < q || q == -1u) {
        printf("incorrect polynomial degrees: %d %d\n", p, q);
        return NT(0);
    }
    NT v1lc = v1[p], invlc = v1lc.pow(zmod::modulus - 2);
    // provided that n >= m
    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));

    // [a, b] = G; [c, d] = H
    for(i = 0; i <= p; i++) {
        a[i] = v1[p - i] * invlc;

        if(i <= q) { // p >= q
            b[i] = v2[q - i];
        }   
    }
    d[q] = NT(1); c[0] = NT(1);

    NT res(1), denom(1);
    NT la(1), lc(1);
    j = 0; 
#if 1
      for(j = 0; j < q; j++) {

        NT b0 = b[j];
        for(i = j + 1; i < p + j + 1; i++) { // update p elements
            b[i] = b[i] - a[i] * b0;
        }
        res = res * v1lc;

        // NOTE: not all coefficients need to be shifted down in fact..
        for(i = r - 2; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }

        c[2*q - j] = b0; // no need to shift 'c' explicitly since it only
                         // stores b[j]'s
        a[j] = 0; c[j] = 0;
        c[q] = 0; // with this we emulate double down-shift
    }
#endif
   
// printf("\n++++++++++++++++++ beginning of main iters G:\n");
//         print_vector(a); print_vector(b);
// 
//         printf("\n B:\n");
//         print_vector(c); print_vector(d);
// 
//         printf("\nla = %d lc = %d\n", la.x, lc.x);
//         printf("\n=========================\n");

    for(; j < r; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];
        for(i = j; i < r; i++) {

            NT ai = la*(a[i]*c0 + b[i]*d0),
               bi = lc*(a[i]*b0 - b[i]*a0);
            a[i] = ai, b[i] = bi;
        }

        // r-th entries of c & d are always zeros
        for(i = j; i < r; i++) {
            NT ci = lc*(c[i]*a0 + d[i]*b0),
               di = la*(c[i]*d0 - d[i]*c0);
            c[i] = ci, d[i] = di;
        }

        lc = lc*lc*la; la = a[j];
        res=res*c[j], denom = denom*lc;

        if(denom.x == 0) { // TODO: early quit ??
           printf("ill-conditioned matrix at j = %d\n", j);
            failed = true;
            return NT(0);
        }

        for(i = r - 2; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }
        // the last entry of 'c' is not shifted down
        for(i = r - 2; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[j] = 0; c[j] = 0; //c[q] = 0; // no need to zero explicitly
    }

    NT inv(mod_inverse(denom.x, zmod::modulus));
    return NT(res * inv);
}

template < class NT >
void sylvester_QR_gcd_direct(const std::vector< NT >& v1,
        const std::vector< NT >& v2, unsigned *out = 0,
            unsigned o_stride = 1) {

    // polynomial degrees
    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;

    if(n < m) {
        printf("sylvester_QR_gcd: incorrect parameters !!\n");
        throw 1;
    }
//     printf("n = %d; m = %d\n", n, m);
    NT v1lc = v1[n], invlc(1);//v1lc.pow(zmod::modulus - 2);

    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0)),
            L(r, NT(0));
    for(i = 0; i <= n; i++) {
        a[i] = (i <= n ? v1[n - i]*invlc : NT(0));
        b[i] = (i <= m ? v2[m - i] : NT(0));
        
        if(i < n)
            c[i + m] = v1[n - i]*invlc;
        if(i < m)
            d[i + n] = v2[m - i];
    }

    NT res(1), denom(1);
#if 1
    for(j = 0; j < m; j++) { // first m steps - work with 2 columns only

//         printf("######### after %d iteration(s) #########\n", j);
//         print_vector(a);
//         print_vector(b);
//         print_vector(c);
//         print_vector(d);

        // extract the current top row: g = G[j]
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

        // compute the L-factor: G.J.g^+
        for(i = j; i < r - (m-1) + j; i++) { // update only n + 1 elements
            L[i] = a0 * a[i] + b0 * b[i];
        }

        // when the rows become linearly-dependent: L-factor is the gcd
//         printf("iteration: %d;\n", j);
//         print_inv_vector(L, L[j], j, r);
        NT det = L[j]; // new denominator
       
        // compute (Z - In).L (emulate down-shift)
        for(i = r - 1; (int)i >= j + 1; i--) {
            L[i] = L[i] - L[i - 1];
        }
        // in fact start update from j+1 because a[j] = b[j] = 0
        unsigned k = r - (m-1) + j+1;
        for(i = j; i < std::min(k, r); i++) { 
            a[i] = a[i] * det - L[i] * a0;
            b[i] = b[i] * det - L[i] * b0;
            // can use another thread to compute la = la * det implicitly
        }               

        printf("\n******** j = %d\n", j);
        print_inv_vector(L, L[j], j, r);

        denom = denom * res * res;
        // do not update c & d: just collect the determinant
        res = res * det;
    }
    for(i = j; i < r; i++) { // premultiply columns c & d by det
        c[i] = c[i] * res;
        d[i] = d[i] * res;
    }
#else
    j = 0;
#endif

    for(; j < r; j++) { // next n steps are full steps

        // extract the current top row: g = G[j]
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

//         printf("iteration: %d;\n", j);

        // compute the L-factor: G.J.g^+
        for(i = j; i < r; i++) {
            L[i] = a0 * a[i] + b0 * b[i] - c0 * c[i] - d0 * d[i];
        }

        // when the rows become linearly-dependent: L-factor is the gcd
//         print_inv_vector(L, L[j], j, r);
        NT det = L[j]; // new denominator
    
        denom = denom * res * res;
        res = res * det;

        std::vector< NT > L2(L);     
        // compute (In - Z).L (emulate down-shift)
        for(i = r - 1; (int)i >= j + 1; i--) {
            L2[i] = (L2[i] - L2[i - 1]);
        }

        // compute a[j + 1] and premultiply
        for(i = j; i < r; i++) { // finally compute: G - (In - Z).L.g / det
            a[i] = a[i] * det - L2[i] * a0;
            b[i] = b[i] * det - L2[i] * b0;
            c[i] = c[i] * det - L2[i] * c0;
            d[i] = d[i] * det - L2[i] * d0;
        }
        printf("\n******** j = %d\n", j);
        print_inv_vector(L, L[j], j, r);


        if(a[j+1]==c[j+1] && b[j+1]==d[j+1]) {
            printf("\n************************* gcd computed\n");
            print_inv_vector(L, L[j], j, r);
//             print_inv_vector(L, NT(1), j, r);
            return;
        }
    }
    printf("res = %d; denom = %d\n", res.x, denom.x);
}

template < class NT >
NT sylvester_gcd(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed,
            const limb *testit = NULL) {

    // polynomial degrees
    unsigned p = v1.size()-1, q = v2.size()-1, r = p + q, i, j;
    if(p < q || q == -1u) {
        printf("incorrect polynomial degrees: %d %d\n", p, q);
        return NT(0);
    }
    unsigned m = zmod::modulus;
    NT v1lc = v1[p], invlc = v1lc.pow(zmod::modulus - 2);
    // provided that n >= m
    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));
    std::vector< NT > aa(r), cc(r);

    // [a, b] = G; [c, d] = H
    for(i = 0; i <= p; i++) {
        a[i] = v1[p - i] * invlc;

        if(i <= q) { // p >= q
            b[i] = v2[q - i];
        }   
    }
    d[q] = NT(1); c[0] = NT(1);

    printf("p: %d; q: %d\n", p, q);

      printf("\n================== original: G:\n");
      print_vector(a); print_vector(b);

      printf("\n B:\n");
      print_vector(c); print_vector(d);


    NT res(1), denom(1);
    NT la(1), lc(1);
    j = 0;    
#if 1
    for(j = 0; j < q; j++) {

        NT b0 = b[j];
        for(i = j + 1; i < p + j + 1; i++) { // update p elements
            b[i] = b[i] - a[i] * b0;
        }
        res = res * v1lc;

        // NOTE: not all coefficients need to be shifted down in fact..
        for(i = r - 2; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }

        c[2*q - j] = b0; // no need to shift 'c' explicitly since it only
                         // stores b[j]'s
        a[j] = 0; c[j] = 0;
        c[q] = 0; // with this we emulate double down-shift
    }
      printf("\n================== beginning iter: G:\n");
      print_vector(a); print_vector(b);

      printf("\n B:\n");
      print_vector(c); print_vector(d);
#endif

    printf("a[j]=%d c[j]=%d d[j]=%d\n", a[j].x, c[j].x, d[j].x);

    for(; j < r; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];
        
        for(i = j; i < r; i++) {

            NT ai = la*(a[i]*c0 + b[i]*d0),
               bi = lc*(a[i]*b0 - b[i]*a0);
            a[i] = ai, b[i] = bi;
        }

        // r-th entries of c & d are always zeros
        for(i = j; i < r; i++) {
            NT ci = lc*(c[i]*a0 + d[i]*b0),
               di = la*(c[i]*d0 - d[i]*c0);
            c[i] = ci, d[i] = di;
        }

        printf("\n================== %d iter: G:\n", j);
        print_vector(a); print_vector(b);

        printf("\n B:\n");
        print_vector(c); print_vector(d);

        printf("\nla = %d lc = %d\n", la.x, lc.x);
        printf("\n=========================\n");

        lc = lc*lc*la; la = a[j];
        res=res*c[j], denom = denom*lc;

        a0 = a[j], c0 = c[j], b0 = b[j+1], d0 = d[j+1];
//         printf("%x %x %x %x\n", a0.x,c0.x,b0.x,d0.x);   
        
        NT prod = la*(a0*c0 + b0*d0);    
        printf("prod2: %x\n", prod.x);

        if(prod.x == 0)
             break;

        NT inva(mod_inverse(la.x, m)), invc(mod_inverse(lc.x, m));
        for(i = j; i < r; i++) {
           aa[i] = a[i] * inva, cc[i] = c[i] * invc;
        }
//         print_vector(aa); print_vector(cc);
        
        if(denom.x == 0) { // TODO: early quit ??
           printf("ill-conditioned matrix at j = %d\n", j);
            failed = true;
            return NT(0);
        }

        for(i = r - 2; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }
        // the last entry of 'c' is not shifted down
        for(i = r - 2; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[j] = 0; c[j] = 0; c[q] = 0; // no need to zero explicitly
    }

    // we stop immediately when la == 0 meaning that we hit the gcd !!
    NT inva(mod_inverse(la.x, m)), invc(mod_inverse(lc.x, m));
    for(i = j; i < r; i++) {
       aa[i] = a[i] * inva, cc[i] = c[i] * invc;
    }

//     a0 = a[j] * inva; c0 = c[j] * invc;
    printf("\nstop at %d iter: \n", j);
//     printf("a0: %d c0: %d\n", a0.x, c0.x);
    print_vector(aa, true); print_vector(cc, true);

    NT inv(mod_inverse(denom.x, zmod::modulus));
    res = res * inv;

    return res;
}

} // namespace internal

} // namespace CGAL

#endif // _MATRIX_ALGEBRA_C_
