
//!////////////////////////////// BFI versions ///////////////////////////////

extern unsigned g_max_err;
extern CORE::BigFloat g_smax, g_smin;

template < class NT >
void dump_bfi(const char *name, const NT& x, unsigned i = -1u) {
//     if(i != -1u)
//         std::cout << i << ": ";
//     std::cout << name << ": " << x << "\n";
}

template < >
void dump_bfi< CORE::BigFloat >(const char *name, const CORE::BigFloat& x,
                unsigned i) {
//     if(i != -1u)
//         std::cout << i << ": ";    
//     std::cout << name << ": " << x << "; " << x.err() << "; "
//             << CORE::bitLength(x.m()) << "\n";
//     std::cout << name << ": " << x.m() << "; " << x.err() << "; "
//             << x.exp() << "\n";
    g_max_err = std::max(g_max_err, (unsigned)x.err());
}

template < class NT, class BFI >
NT resultant_modular_sylvester_bfi(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed, BFI) {

    // polynomial degrees
    unsigned p = v1.size()-1, q = v2.size()-1, r = p + q, i, j;
    if(p < q || q == -1u) {
        printf("incorrect polynomial degrees: %d %d\n", p, q);
        return NT(0);
    }
    NT v1lc = v1[p], invlc;
    // = v1lc.pow(zmod::modulus - 2);
    invlc = NT(1) / v1lc;

    std::cout << "v1lc = " << v1lc << "\n";

    // provided that n >= m
    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));

    // [a, b] = G; [c, d] = H
    for(i = 0; i <= p; i++) {
        a[i] = v1[p - i] * invlc;

        if(i <= q) { // p >= q
            b[i] = v2[q - i];
        }   
    }
    d[q] = NT(1);
    
    NT res(1), denom(1);
    NT la(1), lc(1);

#if 1
      for(j = 0; j < q; j++) {

//         printf("\n================== beginning of %d iter: G:\n", j);
//         print_vector(a); print_vector(b);
// 
//         printf("\n B:\n");
//         print_vector(c); print_vector(d);
// 
//         printf("\nla = %d lc = %d\n", la.x, lc.x);
//         printf("\n=========================\n");

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
#else
    for(j = 0; j < q; j++) {
        res = res * v1lc; // NOTE: does it make sense to do this
                          // for the case of monic polynomials ??
    }
#endif
/*
printf("\n++++++++++++++++++ beginning of main iters G:\n");
        print_vector(a); print_vector(b);

        printf("\n B:\n");
        print_vector(c); print_vector(d);

        printf("\nla = %d lc = %d\n", la.x, lc.x);
        printf("\n=========================\n");*/

//     c[0] = 1; j = 0;
    for(; j < r; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

#if 1
        if(!CGAL::zero_in(a0) && !CGAL::zero_in(c0)) {
//         if(!CGAL::is_zero(a0) && !CGAL::is_zero(c0)) {
            NT inva = b0 / a0, invc = d0 / c0;
            NT det =  1 / (1 + inva * invc);
//                  (a0*c0) / (a0*c0 + b0*d0);

//             std::cout << "CASE 1\n";
            for(i = j; i < r; i++) {

                // NOTE: can also alternate divisions here..
                NT ai = (a[i] + b[i]*invc), // divide by det ??
                   bi = (a[i]*inva - b[i]),
                   ci = c[i] + d[i]*inva,
                   di = c[i]*invc - d[i];

                if(j % 2) {
                    ai *= det; bi *= det;
                } else {
                    ci *= det; di *= det;
                }
                a[i] = ai, b[i] = bi, c[i] = ci, d[i] = di;
            }
        } else if(!CGAL::zero_in(c0)) {
//             std::cout << "CASE 2\n";

            NT invd = c0 / d0;
            for(i = j; i < r; i++) {
                NT ai = a[i] * invd + b[i], // divide by det ??
                   bi = a[i],
                   ci = d[i],
                   di = c[i] - d[i] * invd;
                a[i] = ai, b[i] = bi, c[i] = ci, d[i] = di;
            }
        } else { //if(!CGAL::zero_in(a0)) {
//             std::cout << "CASE 3\n";
    
            NT invb = a0 / b0;
            for(i = j; i < r; i++) {
                NT ai = b[i], // divide by det ??
                   bi = a[i] - b[i]*invb,
                   ci = c[i]*invb + d[i],
                   di = c[i];
                a[i] = ai, b[i] = bi, c[i] = ci, d[i] = di;
            }
        }
#else
        
        std::cout << j << " ###################### NEW ITER\n";

        NT det = (a0*c0 + b0*d0); // == a[j]
//         dump_bfi("a0", a0);
//         dump_bfi("b0", b0);
//         dump_bfi("c0", c0);
//         dump_bfi("d0", d0);
//         dump_bfi("det", det, j);

        for(i = j; i < r; i++) {

            NT ai = (a[i]*c0 + b[i]*d0),
               bi = (a[i]*b0 - b[i]*a0),
               ci = (c[i]*a0 + d[i]*b0),
               di = (c[i]*d0 - d[i]*c0);

            if(j % 2 == 0) {
                ai /= det; bi /= det;
//                 ai.makeExact(); bi.makeExact();
            } else {
                ci /= det; di /= det;
//                 ci.makeExact(); di.makeExact();
            }
            a[i] = ai, b[i] = bi, c[i] = ci, d[i] = di;
        }
#endif
        std::cout << "===================================================\n";
        bool inc = false;
        for(i = j; i < r; i++) {
            dump_bfi("in loop a", a[i], i);
            dump_bfi("in loop b", b[i], i);
            dump_bfi("in loop c", c[i], i);
            dump_bfi("in loop d", d[i], i);
        }

        res = res * a[j] * c[j]; // looks like a[j] & det cancel out

        for(i = r - 2; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }
        // the last entry of 'c' is not shifted down
        for(i = r - 2; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[j] = 0; c[j] = 0;
        c[q] = 0; // no need to zero explicitly

//         printf("\n================== %d iter: G:\n", j);
//         print_vector(a); print_vector(b);
// 
//         printf("\n B:\n");
//         print_vector(c); print_vector(d);
//         printf("\n=========================\n");

    }
    
//!////////////////////////////// required while no inverse is computed
//    denom = zmod(1);
//!//////////////////////////////
//     NT inv(mod_inverse(denom.x, zmod::modulus));
    NT inv = NT(1);

    return NT(res * inv);
}

template < class NT >
void expand_bfi(NT& bfi) {
    
    NT x0 = CGAL::median(bfi);
    
    x1 = Float(Integer(c.err()), 0, c.exp());
}

//! computes the determinant of Sylvester matrix defined by \c v1 and \c v2
//! as well as the \c cfs_s and \c cfs_t (coefficients of s(x) and t(x) in
//! \c res(f,g) = f(x)s(x) + g(x)t(x)
//! \c deg(cfs_s) <\ \c deg(v2) and \c deg(cfs_t) <\ \c deg(v1)
//! \c stride_s/t - output data strides
template < class NT >
NT sylvester_embedding_factorize_bfi2(const std::vector< NT >& v1,
        const std::vector< NT >& v2, limb *cfs_s, limb *cfs_t,
            unsigned stride_s, unsigned stride_t) {

    // polynomial degrees
    unsigned p = v1.size()-1, q = v2.size()-1, r = p + q, i, j;
    if(p < q || q == -1u) {
        printf("incorrect polynomial degrees: %d %d\n", p, q);
        return NT(0);
    }

// TODO TODO TODO: the algorithm does not work unless modular arithmetic is
// used: for bigfloats need to use division in each step, otherwise the
// algorithm hangs

    // provided that n >= m
    std::vector< NT > a(2*r, NT(0)), b(2*r, NT(0)),
        c(r+1, NT(0)), d(r+1, NT(0));

    NT v1lc = v1[p], invlc;
    invlc = NT(1) / v1lc;//v1lc.pow(zmod::modulus - 2);

    // G=(a,b) B=(c,d)
    for(i = 0; i < r; i++) {
        a[i] = (i <= p ? v1[p - i] * invlc : NT(0));
        b[i] = (i <= q ? v2[q - i] : NT(0));
    }
    a[r] = NT(1); b[r + q] = NT(1); // r + p or r + q ??
    c[0] = NT(1); d[q] = NT(1);

    // no need for inverting the leading coefficients.. so simple in fact..
    NT res(1), denom(1);
    NT la(1), lc(1);

      for(j = 0; j < q; j++) {
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

        // b[j] == 0
        for(i = j + 1; i < p + j + 1; i++) {
            b[i] = b[i] - a[i]*b0;
        }
        b[r + j] = -b0;
        c[q] = b0;

        res = res * v1lc;

        //TODO TODO: can avoid explicit shiftings using proper indexing
        for(i = p + j; (int)i >= (int)j; i--) {
            a[i + 1] = a[i]; // shiftdown
        }

        // the last entry of 'c' is not shifted down
        for(i = q + j; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[r] = 0, a[j] = 0;
        c[j] = 0; c[q] = 0; c[r] = 0;

//         printf("\n================== %d iter: G:\n", j);
//         print_vector(a); print_vector(b);
// 
//         printf("\n B:\n");
//         print_vector(c); print_vector(d);
// 
//         printf("\nla = %d lc = %d\n", la.x, lc.x);
//         printf("\n=========================\n");
    }

    a[r + q] = 1; // this is to shift the one

//      printf("++++++++++++++++++ main iterations staring at: %d ++++++++\n", j);

    for(; j < r; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

        NT det = (a0*c0 + b0*d0); // == a[j]

        for(i = j; i < 2*r; i++) {

            NT ai = (a[i]*c0 + b[i]*d0),
               bi = (a[i]*b0 - b[i]*a0);
            if(j % 2 == 0) {
                ai /= det; bi /= det;
            }
            a[i] = ai, b[i] = bi;
        } 

        for(i = j; i < r; i++) {
            NT ci = (c[i]*a0 + d[i]*b0),
               di = (c[i]*d0 - d[i]*c0);

            if(j % 2 == 1) {
                ci /= det; di /= det;
            }
            c[i] = ci, d[i] = di;
        }

        res = res * a[j]*c[j];
#if 0
        // b[r+j] == 0
        for(i = j; i < r + 1 + j; i++) {

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
#endif

        // remark that 'a' is not changes at all - just shifted down..
        // supposedly for the first 'q' iterations
    
        a0 = a[j], c0 = c[j]; // leading elements

        if(j == r - 1) {
            NT inv = a0;
            std::cout << "inv = " << inv << "\ndenom = " << denom << "\n";
            inv = NT(1) / inv;
//             res = res * denom;
            inv = inv * res; // why divide by a0 ??

//             printf("%d: ================ L factor:\n", j);
            limb *pout = cfs_s;
            for(i = r; i < r + q; i++, pout += stride_s) {
        // first q factors multiplied by v1lc to compensate for monic poly
                NT s = a[i] * inv * invlc;
//                 pout[0] = s.x;
                //printf("%d: %d\n", (i-r), s.x);
                std::cout << (i - r) << "; " << s << "\n";
            }
    
            for(pout = cfs_t; i < 2*r; i++, pout += stride_t) {
                NT s = a[i] * inv;
//                 pout[0] = s.x;
//                 if(s.x != invtbl[i-r]) {
//                     printf("WRONG entry (res, truth): %d %d \n",
//                         s.x, invtbl[i-r]);
//                     exit(1);
//                 }
//                printf("%d: %d\n", (i-r), s.x);
                std::cout << (i - r) << "; " << s << "\n";
            }
            break;
        }

        for(i = j + r; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }
        // the last entry of 'c' is not shifted down
        for(i = r - 2; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[r] = 0, a[j] = 0;
        c[j] = 0; c[q] = 0; c[r] = 0;
    }
    return res;
}
//!///////////////////////////////////////////////////////////////////////////

//! computes the determinant of Sylvester matrix defined by \c v1 and \c v2
//! as well as the \c cfs_s and \c cfs_t (coefficients of s(x) and t(x) in
//! \c res(f,g) = f(x)s(x) + g(x)t(x)
//! \c deg(cfs_s) <\ \c deg(v2) and \c deg(cfs_t) <\ \c deg(v1)
//! \c stride_s/t - output data strides
template < class NT >
NT sylvester_embedding_factorize(const std::vector< NT >& v1,
        const std::vector< NT >& v2, limb *cfs_s, limb *cfs_t,
            unsigned stride_s, unsigned stride_t) {

    // polynomial degrees
    unsigned p = v1.size()-1, q = v2.size()-1, r = p + q, i, j;
    if(p < q || q == -1u) {
        printf("incorrect polynomial degrees: %d %d\n", p, q);
        return NT(0);
    }

    // provided that n >= m
    std::vector< NT > a(2*r, NT(0)), b(2*r, NT(0)),
        c(r+1, NT(0)), d(r+1, NT(0));

    NT v1lc = v1[p], invlc = v1lc.pow(zmod::modulus - 2);

    // G=(a,b) B=(c,d)
    NT inv1(1), inv2(1);
    for(i = 0; i < r; i++) {
        a[i] = (i <= p ? v1[p - i] * invlc : NT(0));
        b[i] = (i <= q ? v2[q - i] * inv2 : NT(0));
    }
    a[r] = NT(1); b[r + q] = NT(1); // r + p or r + q ??
    c[0] = NT(1); d[q] = NT(1);

    // no need for inverting the leading coefficients.. so simple in fact..

    // inverse table for debugging
    limb res_truth = 9043741;
    limb invtbl[] = {
        897686, 2227936, 3931761, 15631043, 873029, 13626037, 12730690,
        1434180, 13479767, 15044170, 14720414, 7535206, 226005, 10515173};

    NT res(1), denom(1);
    NT la(1), lc(1);

      for(j = 0; j < q; j++) {
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

        // b[j] == 0
        for(i = j + 1; i < p + j + 1; i++) {
//         for(i = j; i < 2*r; i++) {
            b[i] = b[i] - a[i]*b0;
        }
        b[r + j] = -b0;
        c[q] = b0;

        res = res * v1lc;

        for(i = p + j; (int)i >= (int)j; i--) {
            a[i + 1] = a[i]; // shiftdown
        }

        // the last entry of 'c' is not shifted down
        for(i = q + j; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[r] = 0, a[j] = 0;
        c[j] = 0; c[q] = 0; c[r] = 0;

//         printf("\n================== %d iter: G:\n", j);
//         print_vector(a); print_vector(b);
// 
//         printf("\n B:\n");
//         print_vector(c); print_vector(d);
// 
//         printf("\nla = %d lc = %d\n", la.x, lc.x);
//         printf("\n=========================\n");
    }

    a[r + q] = 1; // this is to shift the one

//     printf("++++++++++++++++++ main iterations staring at: %d ++++++++\n", j);

    for(; j < r; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

//         printf("a[r+j]=%d b[r+j]=%d c[r-1]=%d d[r-1]=%d",
//             a[r+j].x, b[r+j].x, c[r-1].x, d[r-1].x);

        // b[r+j] == 0
        for(i = j; i < r + 1 + j; i++) {

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

        // remark that 'a' is not changes at all - just shifted down..
        // supposedly for the first 'q' iterations
    
        a0 = a[j], c0 = c[j]; // leading elements

        if(j == r - 1) {
            NT inv = a0;
            inv = inv.pow(zmod::modulus - 2);
            denom = denom.pow(zmod::modulus - 2);
            res = res * denom;

            inv = inv * res;

//             printf("%d: ================ L factor:\n", j);
            limb *pout = cfs_s;
            for(i = r; i < r + q; i++, pout += stride_s) {
        // first q factors multiplied by v1lc to compensate for monic poly
                NT s = a[i] * inv * invlc;
                pout[0] = s.x;
            }
    
            for(pout = cfs_t; i < 2*r; i++, pout += stride_t) {
                NT s = a[i] * inv;
                pout[0] = s.x;
            
//                 if(s.x != invtbl[i-r]) {
//                     printf("WRONG entry (res, truth): %d %d \n",
//                         s.x, invtbl[i-r]);
//                     exit(1);
//                 }
//                  printf("%d ", s.x);
            }
            break;
        }

        for(i = j + r; (int)i >= (int)j; i--) { // shiftdown of r + 1 elems
            a[i + 1] = a[i]; // shiftdown
        }
        // the last entry of 'c' is not shifted down
        for(i = r - 2; (int)i >= (int)j; i--) {
            c[i + 1] = c[i]; // shiftdown
        }
        // with this we emulate double down-shift
        a[r] = 0, a[j] = 0;
        c[j] = 0; c[q] = 0; c[r] = 0;

    }
//     printf("\n\nresultant: %d\n", res.x);
//     if(res.x != res_truth)
//         printf("WRONG resultant !!\n");

    return res;
}

template < class NT >
void update_vectors(const std::vector< NT >& xs, std::vector< NT >& a,
         unsigned j, unsigned n) {

    unsigned i;
        NT x0 = xs[j];
        for(i = j; i < n; i++) {
            a[i] = (xs[i] - x0) * a[i];
        }

// after iteration j last non-zero element is a[n+j]
        for(i = n*2 - 1; i >= n + 2; i--) { // n-1 times
            a[i] = a[i - 2] - x0 * a[i]; // to protect from rewriting
        }

        a[n] = NT(0) - x0 * a[n];
        a[n + 1] = NT(0) - x0 * a[n + 1];
}

//! computes resultant cofactors ???
//! \c r - resultant of f and g
//! \c fs - values of f at points \c xs , \c gs - values of g at points \c xs
//! \c n - number of interpolation points
template < class NT >
void quasi_vandermonde_interpolate(const std::vector< NT >& xs,
        const std::vector< NT >& fs, const std::vector< NT >& gs, NT r,
        unsigned n) {

    std::vector< NT > a(2*n, NT(0)), b(2*n, NT(0)), c(2*n, NT(0));

    // G = (a, b, c) B = (d, e, f)
    printf("-----inputs (xs, fs, gs):\n");
    print_vector(xs, true);
    print_vector(fs, true);
    print_vector(gs, true);

    for(unsigned i = 0; i < n; i++) {
        a[i] = fs[i], b[i] = gs[i];
    }
    a[n] = NT(1), b[n + 1] = NT(1);

    printf("\n-----------------original G:\n");
    print_vector(a); print_vector(b); print_vector(c);
    printf("---------------\n\n");

    unsigned i, j;
    NT ld(1);

    NT a0 = a[0], b0 = b[0], c0;
    j = 0;

    // you know that c[n] is trivially one
    for(i = 1; i < n; i++) { // updating [a,b,c] (generator G)
        NT bi = a0 - a[i], ci = b[i]*a0 - a[i]*b0;
        b[i] = bi, c[i] = ci;
    }
    b[n] = -NT(1), c[n] = -b0; // n+j
    b[n+1] = 0, c[n+1] = a0; // n+j+1

//     for(i = j; i < n + 1; i++) {
//         NT di = (d[i]*a0 + e[i]*b0 + f[i]*c0),
//            ei = f[i], fi = e[i];
//         d[i] = di, e[i] = ei, f[i] = fi;
//     }

    printf("\nafter %d th iteration: G: (before update)\n ", j);
        print_vector(a); 

    update_vectors(xs, a, j, n);

    printf("\nafter %d th iteration: G (after update):\n", j);
    print_vector(a); print_vector(b); print_vector(c);
    printf("---------------\n\n");

    //NOTE: if a0 == 0 algorithm breaks down..
    ld = a0;
    NT e0(0), f0(1);

    for(j = 1; j < n; j++) {
 
        NT a0 = a[j], b0 = b[j], c0 = c[j];

//         printf("a[n+j]=%d b[n+j]=%d c[n+j]=%d\n",
//             a[n+j].x, b[n+j].x, c[n+j].x);
//         printf("a[n+j+1]=%d b[n+j+1]=%d c[n+j+1]=%d\n",
//             a[n+j+1].x, b[n+j+1].x, c[n+j+1].x);

        for(i = j + 1; i < n + j; i++) { // updating [a,b,c] (generator G)
//! here the trick is that e0 initially is zero, then it becomes non-zero
            NT ai = (b[i]*e0 + c[i]*f0),
               bi = (b[i]*c0 - c[i]*b0), ci = (b[i]*a0 - a[i]*b0);
            a[i] = ai, b[i] = bi, c[i] = ci;
        }
        NT ai = c[n+j]*f0;
        b[n+j] = -c[n+j]*b0;
        c[n+j] = -a[n+j]*b0;
        a[n+j]=ai;

        // NOTE: you can group and compute c[n+j+1] and a[j]
        b[n+j+1] = 0;
        c[n+j+1] = -a[n+j+1]*b0;
        a[n+j+1]=0;

        a[j] = (b[j]*e0 + c[j]*f0);
        //! NOTE: you can multiply e0 & f0 by any factor

//             e[j+1] = mul*f0*a0, f[j+1] = NT(0) - mul*(b0*e0+c0*f0);
        // e0 & f0 somehow play a role of "remembering" old values..
            e0 = f0*a0; f0 = -a[j];

//         for(i = j+2; i < n+1; i++) {//n + 1; i++) {
//             NT di = (d[i]*a0 + e[i]*b0 + f[i]*c0),
//                 ei = f0*(d[i]*a0 + e[i]*b0) - f[i]*b0*e0,
//                 fi = NT(0) - d[i]*(b0*e0+c0*f0);
//             d[i] = di, e[i] = ei, f[i] = fi;
//         }
//         printf("%d-th: before updating ld: factor=%d rev=%d, ld/a0=%d\n",
//                 j, factor.x, (NT(0)-factor).x, (ld*invld).x);
        if(j == n-1)
            ld = ld*a[j];

        printf("\nafter %d th iteration: G: (before update)\n ", j);
        print_vector(a); 

        update_vectors(xs, a, j, n);

        printf("\nafter %d th iteration: G:\n", j);
        print_vector(a); print_vector(b); print_vector(c);

        printf("\nld = %d\n", ld.x);
        printf("---------------\n\n");
    }

    NT inv(ld);
    inv = inv.pow(zmod::modulus - 2);
    e0 = -r;

    printf("e0=%d ld=%d\n", e0.x, ld.x);

    Simple_vector< NT > res(n);
    for(j = 0; j < n; j++) {
        // in fact d0 == f0 == 0
        res[j] = (b[j + n]*e0)*inv;
    }
    printf("final solution:\n");
    print_vector(res);

    for(j = 0; j < n; j++) {
        // in fact d0 == f0 == 0
        res[j] = NT(0)-res[j];
    }
    printf("reversed solution:\n");
    print_vector(res);

}

template < class NT >
NT sylvester_QR_gcd_givens(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed) {

    // polynomial degrees
    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;

    if(n < m) {
        printf("sylvester_QR_gcd: incorrect parameters !!\n");
        throw 1;
    }

    printf("n = %d; m = %d\n", n, m);
    NT v1lc = v1[n], invlc(1);//v1lc.pow(zmod::modulus - 2);

    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(r, NT(0)), d(r, NT(0));
    for(i = 0; i <= n; i++) {
        a[i] = (i <= n ? v1[n - i]*invlc : NT(0));
        b[i] = (i <= m ? v2[m - i] : NT(0));
        
        if(i < n)
            c[i + m] = v1[n - i]*invlc;
        if(i < m)
            d[i + n] = v2[m - i];
    }

    printf("\n## original:\n");
    print_vector(a, true);
    print_vector(b, true);
    print_vector(c, true);
    print_vector(d, true);

    NT det(1), denom(1), la(1);
#if 0
    for(j = 0; j < m; j++) { // first m steps - work with 2 columns only
                    
        NT a0 = a[j], b0 = b[j];
        for(i = j; i < r - (m-1) + j; i++) { // only n+1 element needs to be
                                             // updated
            NT ai = a0 * a[i] + b0 * b[i];
            b[i] = (a0 * b[i] - b0 * a[i]);
            a[i] = ai;
        }
        la = la * a[j];
        printf("######## iter: %d\n", j);
        print_inv_vector(a, a[j], j, r);

        for(i = r - 2; (int)i >= (int)j; i--) { 
            a[i + 1] = a[i]; // shiftdown
        }
    }
#else
    j = 0;
#endif
    printf("\n# of d_null iterations: %d\n", n-m);

    NT mb(1), mc(la), md(la);   // scaling factors for b & c columns
    for(; j < r; j++) { // next n steps are full steps

        printf("## main iteration: %d\n", i);

        NT a0 = a[j], b0 = b[j], c0 = c[j], d0 = d[j];

        NT mb0 = mb * b0; // eliminate column b
        for(i = j; i < r; i++) {
            NT ai = a0 * a[i] + mb0 * b[i];
            b[i] = a0 * b[i] - b0 * a[i]; 
            a[i] = ai;
        }
        mc = mc * a[j], md = md * a[j]; // at each step scale factors
                            // make the denominators equal for each column
        // use one more thread which has c0 == 0
        // and a[i] == mc, a[j] == md -> to make updates automatic

        if(c0 != NT(0)) {
            a0 = a[j]; // update a0
            NT mc0 = -mc * c0; //! hyperbolic rotation !!! 
            for(i = j; i < r; i++) {
                NT ai = a0 * a[i] + mc0 * c[i];
                c[i] = a0 * c[i] - c0 * a[i];
                a[i] = ai;
            }
            mb = mb * a[j], md = md * a[j];
        } else printf("############ c null ###########\n");

        if(d0 != NT(0)) {
            a0 = a[j];

            NT md0 = -md * d0; //! hyperbolic rotation !!! 
            for(i = j; i < r; i++) {
                NT ai = a0 * a[i] + md0 * d[i];
                d[i] = a0 * d[i] - d0 * a[i];
                a[i] = ai;
            }
            mb = mb * a[j], mc = mc * a[j]; 
        } else printf("############ d null ###########\n");

        print_inv_vector(a, a[j], j, r);

        for(i = r - 2; (int)i >= (int)j; i--) { 
            a[i + 1] = a[i]; // shiftdown
        }
    }

    return NT(1);
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

    NT res(1), denom(1);
    NT la(1), lc(1);
    j = 0;    
#if 0
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


//! annihilates an element \c b[b_i] using columns \c a with leading
//! element \c a[a_i] . The columns  are of length \c n . 
//! New scale factors are returned in \c la and \c lb .
template < class NT >
void givens_rotate(std::vector< NT >& a, std::vector< NT >& b,
        unsigned a_i, unsigned b_i, unsigned n, NT& la, NT& lb,
        bool hyperbolic = false) {

    NT a0 = a[a_i], b0 = b[b_i], lb_m_a0 = lb * a0, la_m_b0 = la * b0;
        
    if(hyperbolic)
        la_m_b0 = NT(0) - la_m_b0;

    for(unsigned j = 0; j < n; j++) {
        NT aj = lb_m_a0 * a[j + a_i] + la_m_b0 * b[j + b_i];
        b[j + b_i] = a0 * b[j + b_i] - b0 * a[j + a_i];
        a[j + a_i] = aj;
    }
    la = la * lb * a[a_i], lb = a[a_i];
}


template < class NT >
NT sylvester_QR_gcd_old(const std::vector< NT >& v1,
        const std::vector< NT >& v2, bool& failed) {

    std::cout << "\n================== sylvester_QR_gcd\n";

    unsigned n = v1.size()-1, m = v2.size()-1, r = n + m, i, j;
    if(n < m || m == -1u) { 
        printf("incorrect polynomial degrees: %d %d\n",n, m);
        return NT(0);
    }
        
    NT v1lc = v1[n], invlc(1);// = v1lc.pow(zmod::modulus - 2);

    std::vector< NT > a(r, NT(0)), b(r, NT(0)), c(n, NT(0)), d(n, NT(0));
    for(i = 0; i < r; i++) {
        a[i] = (i <= n ? v1[n - i]*invlc : NT(0));
        b[i] = (i <= m ? v2[m - i] : NT(0));
    }
    for(i = 0; i < n; i++) {
        c[i] = v1[n - i] * invlc;
        d[i] = (i >= n - m ? v2[n - i] : NT(0));
    }

    NT det(1), denom(1), la(1), lb(1);
    for(i = 0; i < m; i++) { // first m steps - annihilate only 
                            // 2nd column elements 
        
        // instead of downshift we continue using a[0] but
        // shift down b's by means of indexing 
        givens_rotate<NT>(a, b, 0, i, r - i, la, lb);
        
        printf("## iter: %d\n", i);
        print_inv_vector(a, a[0], 0, r - i);
        det = det * a[0], denom = denom * la;
    }

// TODO: rewrite algorithm in a normal way without complicated
// indexing
    NT lc(1), ld(1);
    for(; i < r; i++) { // next n steps are full steps

        unsigned k = r - i, s;

        NT a0 = a[0], b0 = b[i], lb_m_a0 = lb * a0, la_m_b0 = la * b0;
        NT c0 = c[i - m], d0 = d[i - m],
            ld_m_c0 = ld * c0, lc_m_d0 = lc * d0;

//         printf("## main iteration: %d %d %d %d\n", i, a0.x, b0.x,
//                 c0.x, d0.x);

        for(j = 0; j < k; j++) {
            NT t = lb_m_a0 * a[j] + la_m_b0 * b[j + i];
            b[j + i] = a0 * b[j + i] - b0 * a[j], a[j] = t;

            s = j + i - m;
            t = ld_m_c0 * c[s] + lc_m_d0 * d[s];
            d[s] = c0 * d[s] - d0 * c[s], c[s] = t;
        }
        la = la * lb * a[0], lb = a[0];
        lc = lc * ld * c[i - m], ld = c[i - m];

        givens_rotate(a, c, 0, i - m, k, la, lc, true);
        // vectors have the form: [ a / sqrt(la) ] etc..

//         print_vector(a);
//         print_vector(b);
//         print_vector(c);
//         print_vector(d);

//         printf("## main iteration: %d\n", i);
//         print_inv_vector(a, a[0], 0, k);

//         NT inva(mod_inverse(a[0].x, zmod::modulus)),
//            invc(mod_inverse(lc.x, zmod::modulus));
//         std::vector< NT > aa(r, NT(0)), cc(r, NT(0));
// 
//         for(j = 0; j < k; j++) {
//             aa[j] = a[j] * inva;//, cc[j] = c[j] * c[0] * invc;
//         }
//         // a^2 / la to get the value
//         print_vector(aa, true);

        det = det * a[0];
        denom = denom * la;
    }

    //zmod root = denom.pow((zmod::modulus+1)>>2);
    std::cout << det << "/sqrt(" << denom << ")\n\n";

    NT res;
//     mpz_sqrt(res.get_mp(), denom.get_mp());
//     std::cout << "sqrt = " << res << "\n\n";
//     std::cout << (det / res) << "\n\n";

    return NT(1);
}