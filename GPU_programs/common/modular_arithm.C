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

#include <include/modular_arithm.h>

namespace CGAL {

namespace internal {

#if 1

void make_oddprime_bitarray(ulong n, Bit_vector& ba) {

    ba = Bit_vector(n/2, 1); 
    const ulong m = ba.size();

    ba[0] = 0;
    ulong r = (unsigned)sqrt((double)n);
    ulong i = 3, ih = i/2;
    
    while(i <= r) {
        if(ba[ih]) {
            for(ulong kh = (i*i)/2; kh < m; kh += i)  
                ba[kh] = 0;;
        }
   
        while(++ih < m) // find next set bit
            if(ba[ih])
                break;
        i = 2*ih + 1;
    }
}

// Return next prime >= n.
// Return zero if table of primes is not big enough.
ulong next_small_prime(ulong n) {
    if(n <= 2)  
        return 2;

    ulong nh = n/2, m = oddprime_array.size();
    if(nh >= m)
        return 0;

    while(nh < m) {
        if(oddprime_array[nh])
            break;
        nh++;
    }

    if(nh == m)
        return 0;
    return (2 * nh + 1);
}

bool is_prime(umod_t x) {

    if((x & 1) == 0)
        return (x == 2);

    umod_t xh = x / 2;
    if(xh < (umod_t)oddprime_array.size())
        return oddprime_array[xh];

    if(x % 3 == 0 || x % 5 == 0 || x % 7 == 0 || x % 11 == 0 ||
            x % 13 == 0)
        return false;

    unsigned prime = 15;
    while(1) {
        unsigned next = next_small_prime(prime + 2);
        if(next == 0)
            break;
        if((x % next) == 0)
            return false;
        prime = next;
    }
    unsigned sq = (unsigned)floor(sqrt((double)x));

     for(unsigned c = prime + 2; c < sq; c += 2) {
        if(x % c == 0)
            return false;
    }
    return true;
}
#endif

#if 0
// computes factorization of m as a vector of pairs (prime, exp),
// return true if m is composite
bool factorize(umod_t m, Pair_vector& fact) {

    if(is_prime(m)) {
        fact.push_back(Int_pair(m, 1));
        return false;
    }

    umod_t maxv = (umod_t)sqrt((double)m) + 1, v = 1;
    while(1) {

        umod_t next = next_small_prime(v + 1);
        if(next == 0) {
            // find prime >= next
            next = v + 1;
            if((next & 1) == 0)
                next++;
            while(!is_prime(next)) 
                next += 2;
        }
        if(next > maxv) 
            break;        

        umod_t exp = 0;
        while((m % next) == 0) {
            m /= next;
            exp++;
        }
        if(exp) {
            fact.push_back(Int_pair(next, exp));
            maxv = (umod_t)sqrt((double)m) + 1;
        }
        v = next;
    }

    if(m != 1) {// rest is prime factor
        if(is_prime(m))
        fact.push_back(Int_pair(m, 1));
        else
        printf("attention: %d not a prime!\n", m);
    }
    return true;
}

// Return the order of x (mod m).
umod_t order_mod(umod_t x, umod_t m, const Pair_vector& phifact)
{
  
    umod_t h = m - 1;// phifact.product();
    umod_t e = h;

   // printf("computing ord(%u) mod %u ========= \n", x, m);
    for(Pair_vector::const_iterator pit = phifact.begin(); 
        pit != phifact.end(); pit++) {

            umod_t p = pit->first;
            umod_t f = pit->second;

        f = ipow(p, f); // f = p^f
        //printf("p = %d, f = %d\n", p, f);

        e /= f;
        umod_t g1 = pow_mod(x, e, m);  // Cohen
      
        while(g1 != 1) {
            g1 = pow_mod(g1, p, m);
            e *= p;
        }
    }
    return e;
}

// Return element of maximal order modulo m (m must be a prime) 
// phifact must be the factorization of m - 1
umod_t maxorder_element_mod(umod_t m, umod_t k, const Pair_vector& phifact) {
    if(m == 2)  
        return 1;

    umod_t mo = m - 1;
    for(umod_t x = m-1 ; x >=2; x--) {
        umod_t xo = order_mod(x, m, phifact);
        if(mo == xo) {
            return  x;
        }
    }
    return 0;
}
#endif

#if 1

void get_next_mod(unsigned& startm, unsigned nskips) {

    static int ccnt = 0;
    //! # of moduli of the form 2^9*k + 1: 4229
    //! last modulo: 0x2a01; total # of bits: 94979

    ccnt++;
    unsigned m = startm;
    float nbits = 0;

    while(m != 0) {
        if(is_prime(m)) {
            if(nskips == 0)
                break;
//             nbits += log((float)m) / log(2.0f);
            nskips--;
        }
        m--;
    }
    if(m == 0) {
        printf("FATAL: moduli exhausted..\n");
        exit(1);
    }
    //printf("total moduli: %d; nbits: %f\n", ccnt, nbits);
    zmod::modulus = m;
    startm = m;
}

//! generates moduli table of size \c nmods and saves it to the file \c fname
//! reports error in case moduli exhausted
void generate_mod_table(unsigned nmods, const char *fname) {

    unsigned i, m, nskips = 0;
#if CUMP_USE_32BIT_MODULI_SET
    m = 1U << 31;//0xffffffffU;
#else
    m = 0xffffffU;
#endif

    // one *guard* float for negative offsets
    double *bits = new double[nmods + 1], *pbits = bits + 1;
    bits[0] = 0;

    FILE *fp = fopen(fname, "wb");
    fprintf(fp, "%d ", nmods);
    
    for(i = 0; i < nmods; i++, m--) {

        get_next_mod(m, nskips);
        if(m == 0) {
            printf("FATAL: moduli exhausted..\n");
            exit(1);
        }
        MOD_entry ee; ee.m = m;
        pbits[i] = pbits[i-1] + log2((double)m);
        printf("%d: %x bits: %.10f\n", i, m, pbits[i]);

#if CUMP_USE_32BIT_MODULI_SET
// NOTE NOTE: this generates other different invm!!
        ee.invm_d = (double)(1U << 30) / m;

        uint64 u1, u2, beta = 1ULL << 32;
        egcd(beta, (uint64)m, u1, u2);
        if((int64)u2 < 0) 
            u2 += beta;
        ee.mu = (unsigned)(beta - u2);
#else
        ee.invk = 65536.0f / m, ee.invm = 1.0f / m;

        unsigned beta = 1 << 24, u1, u2;
        egcd(beta, m, u1, u2);
        if((int)u2 < 0) // compute m^-1 mod beta needed for Montgomery inv
            u2 += beta;
        ee.mu = beta - u2;
#endif // CUMP_USE_32BIT_MODULI_SET

        fwrite(&ee, sizeof(MOD_entry), 1, fp);
    }
    fwrite(pbits, sizeof(double), nmods, fp);
    fclose(fp);
    delete []bits;
}

#endif

//! computes x ^ -1 mod m using Montgomery inverse algorithm
//! \c mu = -m ^ -1 mod beta (beta = 2^24 or 2^32 depending on mod bitlength)
unsigned mod_inverse_montgomery(unsigned x, unsigned m, unsigned mu) {

    int Luv = x, Ruv = m, Lrs = 1, Rrs = 0, k = 0;
    while(1) {

        int tmprs = Rrs;
        if(Luv & 1) {
            int safeuv = Luv;

            if((safeuv ^ Ruv) < 0)
                Luv += Ruv;
            else
                Luv -= Ruv;

            if((Luv ^ safeuv) < 0) {
                Ruv = safeuv, tmprs = Lrs;
            }
            Lrs += Rrs;
        }
        Luv >>= 1, k++;
        Rrs = tmprs * 2;

        if(Luv == 0)
            break;
    }
    Lrs = m - Rrs;
    if(Lrs < 0)
        Lrs += m;

    int r = Lrs;

#if CUMP_USE_32BIT_MODULI_SET
    if(k >= 32) { // r = MonPro(r, 1)

        unsigned C = r, Q = C * mu;
        uint64 prod = (uint64)Q * (uint64)m;
        unsigned lo = (unsigned)prod, hi = (unsigned)(prod >> 32);
        lo = lo + C;
        hi += (lo < C); // add carry out: can use addc here
        r = hi, k -= 32;
    }
    if(k == 0)
        return r;

    // r = MonPro(r, 1 << (32 - k))
    unsigned c_small = r << (32 - k), Q = c_small * mu;
    uint64 prod = (uint64)Q * (uint64)m;
    unsigned lo1 = (unsigned)prod, hi1 = (unsigned)(prod >> 32);

    unsigned lo2 = c_small, hi2 = r >> k; // == c_small ??
    lo1 += lo2;
    hi1 += hi2 + (lo1 < lo2);   // can use addc here
    r = hi1;
    
#else
    if(k > 24) {
        unsigned C = r, Q = my_umul24(C, mu), beta = 1 << 24;
        
        unsigned hi, lo;
        mul_24_wide_host(hi, lo, Q, m);
        lo = my_umul24(lo, 1) + C;
        hi = (hi >> 8) + (lo >= beta);
        r = hi;
        k -= 24;
    }

    unsigned c_small = r << (24 - k), Q = my_umul24(c_small, mu);

    unsigned hi1, lo1, hi2, lo2;
    mul_24_wide_host(hi1, lo1, Q, m);

    lo2 = c_small, hi2 = r >> k; // == c_small ??

    lo1 &= 0xffffff;
    unsigned sum = lo1 + (lo2 & 0xffffff);
    hi1 = (hi1 >> 8) + hi2 + (sum >> 24);
    r = hi1;
#endif // CUMP_USE_32BIT_MODULI_SET

    if(r < 0 || r >= m) {
        printf("%x: wrong r: %x\n", m, r); exit(1); }
    return r;
}

unsigned load_moduli_set(std::vector< MOD_entry >& mod_table,
         std::vector< double >& mod_bits, const char *fname) {

    FILE *fp = fopen(fname, "rb");
    if(fp == NULL) {
        printf("FATAL: could not open moduli set: %s\n", fname);
        return -1u;
    }

    unsigned nmods;
    fscanf(fp, "%d ", &nmods);
    if(nmods == 0 || nmods > 50000) 
        goto Lerr_exit;
    
    mod_table.resize(nmods);
    mod_bits.resize(nmods);
    
    if(fread(mod_table.data(), sizeof(MOD_entry), nmods, fp) < nmods) 
        goto Lerr_exit;
        
    if(fread(mod_bits.data(), sizeof(double), nmods, fp) < nmods)
        goto Lerr_exit;

    printf("Loaded moduli set of size: %d\n", nmods);
    fclose(fp);
    return nmods;

Lerr_exit:
    printf("FATAL: corrupt moduli set: %s\n", fname);
    fclose(fp);
    return -1u;
}

//! returns # of moduli required to recover \c nbits
//! returns 0 if moduli set exhausted
unsigned bits_to_moduli(unsigned nbits,
        const std::vector< double >& mod_bits) {

    int li = 0, ri = (int)mod_bits.size() - 1, mi;
    while(li + 1 < ri) {
        mi = (li + ri) / 2;
//         printf("%d: bits: %.13f\n", mi, mod_bits[mi]);
        if(nbits <= mod_bits[mi]) {
            ri = mi;
        } else
            li = mi;
    }
    if(nbits < mod_bits[0])
        return 1;

    if(nbits > mod_bits[ri]) // moduli exhausted
        return 0;
    return (unsigned)ri+1;
}

//! takes a sequence of moduli \c mods and computes decreasing modular
//! inverses, i.e., for n_mods = 4: m4^-1 mod m3, (m4*m3)^-1 mod m2,
//! (m4*m3*m2)^-1 mod m1
//! \c r must have enough space for \c n_mods-1 residues
//! \c stride - memory stride to access moduli
//! precondition: m[0] > m[1] > m[2] > .. > m[n-1]
void mod_inverse_seq(unsigned *r, const unsigned *mods, unsigned n,
        unsigned stride) {

    const unsigned *pm_i = mods + (n - 1) * stride;
    unsigned *m = new unsigned[n];

    unsigned i, j;
    m[n - 1] = pm_i[0], pm_i -= stride;
    for(j = n - 2; (int)j >= 0; j--, pm_i -= stride) {
        r[j] = m[n - 1]; m[j] = pm_i[0];
//         r2[j - 1] = m[0]; // 1 modulus is not needed but it does not harm
    }
    // moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    // inv_mod[n - 1] = 1
    // inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..

    // 1st row: n-2 elements c[i] * m[0]
    // 2nd row: n-3 elements c[i] * m[0] * m[1] ...
    // 3rd row: n-4 elements c[i] * m[0] * m[1] * m[2] ...
//     pr2 = r2 + n - 2;

    for(i = n - 3; (int)i >= 0; i--) {
        for(j = i; (int)j >= 0; j--) {
            zmod::set_mod(m[j]);
            r[j] = mul_mod(r[j], m[i + 1], m[j]);
        }
    }
    // now r[i] = mods[0]*mods[1]*..*mods[i] mod mods[i+1]
    for(unsigned j = 0; j < n - 1; j++) {
        zmod::set_mod(m[j]);
        r[j] = mod_inverse(r[j], m[j]);
    }
    delete []m;
}

//! adds residue \c x modulo \c m to final result \c r
//! \c nr & \c nprod - length of r and prod resp. (in words)
//! \c inv_prod - \c prod^-1 mod m
//! attention: \c r & \c prod must have enough for residue and modulo resp.
void incremental_CRT(limb *r, unsigned& nr, limb *prod, unsigned& nprod,
        limb x, limb m, limb inv_prod, unsigned n_max) {

//TODO: make sure u compile with 32-bit gmp

    const unsigned N_BITS = sizeof(limb) * 8;

    mpz_t mp_r, mp_prod;
    mpz_init2(mp_r, nr * N_BITS);
    memcpy(mp_r->_mp_d, r, nr * sizeof(limb));
    mp_r->_mp_size = nr;

    mpz_init2(mp_prod, nprod * N_BITS);
    memcpy(mp_prod->_mp_d, prod, nprod * sizeof(limb));
    mp_prod->_mp_size = nprod;

    // need to compute [prod]^-1 mod m, should this be passed as a parameter ?
    // if anyway we need this for kernel ?

    mpz_t mp_inv, mp_m, mp_x, t;
    mpz_init(mp_inv);
    mpz_init_set_ui(mp_m, m);
    //mpz_invert(mp_inv, mp_prod, mp_m); // prod**(-1) mod m2
    mpz_init_set_ui(mp_inv, inv_prod);
    mpz_init_set_ui(mp_x, x);
    mpz_init(t);

//     printf("prod: %x; inv_prod: %x m: %x\n", prod[0], inv_prod, m);

    mpz_mul(t, mp_inv, mp_prod); // mpz_mul_ui ??
    mpz_mod_ui(t, t, m);
    if(t->_mp_d[0] != 1) {
        printf("---------------------wrong inverse!!\n");
        exit(1);
    }

    mpz_sub(t, mp_x, mp_r);
    mpz_mod_ui(t, t, m); // x2 - x1 mod m2
    mpz_mul(t, t, mp_inv);
    mpz_mod_ui(t, t, m); // (x2 - x1)*c mod m2
    mpz_mul(t, t, mp_prod); // s * m1
    mpz_add(mp_r, mp_r, t); // x1 + s*m1

    nr = mp_r->_mp_size;
    if(nr > n_max) {
        printf("r size is too big: %d %d\n", nr, n_max);
        exit(1);
    } else 
        memcpy(r, mp_r->_mp_d, nr * sizeof(limb));

    mpz_mul(mp_prod, mp_prod, mp_m); // add m to moduli product
    nprod = mp_prod->_mp_size;
    if(nprod > n_max) {
        printf("prod size is too big: %d %d\n", nprod, n_max);
        exit(1);
    } else
        memcpy(prod, mp_prod->_mp_d, nprod * sizeof(limb));

    mpz_clear(mp_r);
    mpz_clear(mp_inv);
    mpz_clear(mp_m);
    mpz_clear(mp_prod);
    mpz_clear(mp_x);
    mpz_clear(t);
}

/**
//NOTE: we need a prime table to generate moduli
//     printf("building prime table..\n");
//     make_oddprime_bitarray(small_prime_limit, oddprime_array);
    unsigned i, u = (1u << 24), m = u - 2;
    for(i = 0; i < CUMP_N_MODULI; i++, pconst += c_stride, m--) {

#if 1
        get_next_mod(m, 0);
        printf("%d: modulus: %x\n", i, m);
        float invk = 65536.0f / m, invm = 1.0f / m;

        unsigned beta = 1 << 24;
        unsigned u1, u2;
        egcd(beta, m, u1, u2);
        if((int)u2 < 0) // compute m^-1 mod beta needed for Montgomery inv
            u2 += beta;

        pconst[0] = m, pconst[1] = (unsigned&)invk;
        pconst[2] = (unsigned&)invm; pconst[3] = beta - u2;
        Mods[i] = mod_table[i].m; InvKs[i] = (unsigned&)invk;
        Mus[i] = beta - u2;
#else
//         Mods[i] = mod_table[i].m; InvKs[i] = (unsigned&)mod_table[i].invk;
//         Mus[i] = mod_table[i].mu;
#endif
    }
*/

#if 0
void CRA_bases(unsigned *r, unsigned r_stride,  const unsigned *mods,
        unsigned n_mods, unsigned m_stride) {

    mpz_t M, R;
    mpz_init_set_ui(M, mods[0]);
    mpz_init(R);
    unsigned i;

    const unsigned *pm = mods + m_stride;
    for(i = 1; i < n_mods; i++, pm += m_stride) {
        mpz_mul_ui(M, M, pm[0]);
    }

    for(i = 0, pm = mods; i < n_mods; i++, pm += m_stride, r += r_stride) {
    
        mpz_divexact_ui(R, M, pm[0]); // divide out this modulus
        mpz_mod_ui(R, R, pm[0]);
        unsigned ur = mpz_get_ui(R);
        r[0] = mod_inverse(ur, pm[0]);
//         ur = mul_mod(ur, r[0], pm[0]);
//         printf("%d: %x\n", i, ur);
    }

    mpz_clear(R);
    mpz_clear(M);
}

//! computes weights \c w for core function \c CM 
void core_weights(unsigned *w, unsigned w_stride, unsigned CM,
    const unsigned *bases, const unsigned *mods, unsigned n_mods,
    unsigned m_stride) {

    const unsigned *pm = mods + m_stride;
    unsigned i;
    for(i = 0; i < n_mods; i++, pm += m_stride, w += w_stride) {

// NOTE: some weights must be negative to make sure equation about C(M)
// is satisfied: unclear how to choose them..
        unsigned CMw = CM % pm[0]; // make sure this is a proper residue
        w[0] = mul_mod(CMw, bases[i], pm[0]);
//         std::cout << i << ": " << w[0] << "\n";
    }
}

//! computes evaluates core function \c CM for residue set \c r ,
//! CRA \c bases and set of core \c weights 
void core_function(const unsigned *r, unsigned r_stride, unsigned CM,
    const unsigned *weights, const unsigned *bases,
    const unsigned *mods, unsigned n_mods, unsigned m_stride) {

    const unsigned *pm = mods + m_stride;
    unsigned i;

    uint64 cn(0);
    for(i = 0; i < n_mods; i++, pm += m_stride, r += r_stride) {

        uint64 CBi = ((uint64)CM * bases[i] - weights[i]) / pm[0];
        CBi *= r[0];

//         std::cout << i << ": " << cn << "; " << CBi << "\n";
        printf("%d: %#llx - %#llx\n", i, cn, CBi);
        cn += CBi;
        
//         unsigned CMw = CM % pm[0]; // make sure this is a proper residue
//         w[0] = mul_mod(CMw, bases[i], pm[0]);
    }
}

//!\c bases - CRA bases: (M / m[i])^(-1) mod m[i]
void RNS_magnitude(const unsigned *r, unsigned r_stride,
    const unsigned *bases, const unsigned *mods, unsigned n_mods,
        unsigned m_stride) {

    const unsigned *pm = mods + m_stride;
    unsigned i;

//     typedef CORE::BigFloat Float;
    typedef double Float;
//     CGAL::set_precision(Float(), 1000);
    Float dres(0);

    std::cout.precision(25);
    
    for(i = 0, pm = mods; i < n_mods; i++, pm += m_stride, r += r_stride) {

        unsigned x = mul_mod(r[0], bases[i], pm[0]);
//         Rational ss(x, pm[0]);
//         ss -= Rational(12, n_mods);

        // point is when you add another number here
        // the presense of integer part will screw up every'

//         dres += CGAL::convert_to_bfi(ss);
//         dres = dres - dres.intValue(); //CORE::floor(dres);

        dres += (Float)x / (Float)pm[0];// - 12.0 / (Float)n_mods;
        dres = dres - std::floor(dres);

        std::cout << i << ": " << dres << "\n";
    }
//     unsigned x = mul_mod(r[0], bases[n_mods-1], pm[0]);
//     dres = (Float)pm[0] / ((Float)pm[0] * dres + x);
//     std::cout << i << ": " << (1.0/dres) << "\n";
   
}
#endif

void convert_to_RNS(unsigned *r, unsigned r_stride, const mpz_t mp_in,
        const unsigned *mods, unsigned n_mods, unsigned m_stride) {

//     unsigned i;
//     if(

    mpz_t mp_r;
    mpz_init(mp_r);

    const unsigned *pm = mods;
    for(unsigned i = 0; i < n_mods; i++, pm += m_stride, r += r_stride) {
        mpz_mod_ui(mp_r, mp_in, pm[0]);
        r[0] = mpz_get_ui(mp_r);
    }

    mpz_clear(mp_r);
}

//! reduces input \c in by \c n moduli \c mods , writes results to \c r
//! \c r must have enough space for \c n_mods residues
void convert_to_RNS_pure(unsigned *r, unsigned r_stride, const limb *in,
        unsigned n, const unsigned *mods, unsigned n_mods,
                unsigned m_stride) {

    mpz_t mp_in;
    mpz_init2(mp_in, n * sizeof(unsigned) * 8);
    memcpy(mp_in->_mp_d, in, n * sizeof(limb));
    mp_in->_mp_size = n;

    convert_to_RNS(r, r_stride, mp_in, mods, n_mods, m_stride);
    mpz_clear(mp_in);
}

// given n residues modulo mods[i], computes associated mixed radix digits
// c[i] = (m[0]*m[1]*..*m[i-1])^-1 mod m[i]
// assume: mods[0] < mods[1] < .. < mods[n-1]
//! returns in \c y an MRS digit representation of a large integer given by its residues \c x ; \c x and \c y can overlap
void compute_MR_digits(unsigned *y, const unsigned *x, const unsigned *mods,
        const unsigned *c, unsigned n) {

    unsigned i, j;
    unsigned *M = new unsigned[n]; // moduli products

    y[0] = x[0];

    for(i = 1; i < n; i++) {
        zmod::set_mod(mods[i]);
        unsigned t = sub_mod(x[i], x[0], mods[i]);
        y[i] = mul_mod(t, c[i - 1], mods[i]);
        M[i] = //mods[0];
            mul_mod(mods[0], c[i - 1], mods[i]);
    }

    for(i = 2; i < n; i++) {

        for(j = i; j < n; j++) {
            zmod::set_mod(mods[j]);
            unsigned t = mul_mod(y[i - 1], M[j], mods[j]);
            y[j] = sub_mod(y[j], t, mods[j]);
            M[j] = mul_mod(M[j], mods[i - 1], mods[j]);
             //printf("rr: %x\n", mul_mod(M[j], c[i - 1], mods[j]));
        }
    }
    delete []M;
}

// the same as above but everything is stored in reversed order
//! x_stride: in-out data stride
//! m_stride: moduli stride
void compute_MR_digits_rev(unsigned *y, const unsigned *x, unsigned x_stride,
        const unsigned *mods, unsigned m_stride, const unsigned *c,
            unsigned n) {

    unsigned i, j;
    std::vector< unsigned > Mv(n);
    unsigned *M = (unsigned *)Mv.data();

    const unsigned *px = x + (n - 1) * x_stride,
        *pm = mods + (n - 1) * m_stride, *ppx, *ppm;

    unsigned *py = y + (n - 1) * x_stride, *ppy;
    py[0] = px[0];//y[n - 1] = x[n - 1];
    // moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    // inv_mod[n - 1] = 1
    // inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..

    for(i = n - 2, ppx = px - x_stride, ppy = py - x_stride, 
            ppm = pm - m_stride; (int)i >= 0; i--, ppx -= x_stride, 
                    ppy -= x_stride, ppm -= m_stride) {
        zmod::set_mod(ppm[0]);
        unsigned t = sub_mod(ppx[0], px[0], ppm[0]);
        ppy[0] = mul_mod(t, c[i], ppm[0]);
        M[i] = mul_mod(pm[0], c[i], ppm[0]);
    }
//     return;

    //! NOTE: only MR digits from y[n-1..stop_idx] are computed !
    //! the remaining x[stop_idx-1..0] residues are not modified
    //! \c n_ready - # of MR digits being computed
    int stop_idx = 0, n_ready = n - stop_idx;

    py -= x_stride, pm -= m_stride;
    for(i = n - 3; (int)i >= stop_idx; i--, pm -= m_stride, py -= x_stride) {
        for(j = i, ppy = py - x_stride, ppm = pm - m_stride;
                (int)j >= stop_idx; j--, ppy -= x_stride, ppm -= m_stride) {
            zmod::set_mod(ppm[0]);
            unsigned t = mul_mod(py[0], M[j], ppm[0]);
            ppy[0] = sub_mod(ppy[0], t, ppm[0]);
            M[j] = mul_mod(M[j], pm[0], ppm[0]);
        }
    }

#if 0
    pm = mods + (n - 2) * m_stride;
    py = y + (n - 2) * x_stride;
    //! update the remaining \c stop_idx residues using \c n_ready computed
    //! MR digits
    for(j = 0; j < n_ready-1; j++, py -= x_stride, pm -= m_stride) {
        ppm = mods + (stop_idx-1) * m_stride;
        ppy = y + (stop_idx-1) * x_stride;
        for(i = stop_idx-1; (int)i >= 0; i--, ppm -= m_stride,
                    ppy -= x_stride) {
            zmod::set_mod(ppm[0]);
            unsigned t = mul_mod(py[0], M[i], ppm[0]);
            ppy[0] = sub_mod(ppy[0], t, ppm[0]);
            M[i] = mul_mod(M[i], pm[0], ppm[0]);
        }
    }
#endif
/**
    y[n - 1] = x[n - 1];
    // moduli in descreasing order: m[0] > m[1] > ... > m[n-1]
    // inv_mod[n - 1] = 1
    // inv_mod[n - 2] = m[n - 1]^-1 mod m[n - 2]
    // inv_mod[n - 3] = m[n - 1]m[n - 2]^-1 mod m[n - 3] ..

    for(i = n - 2; (int)i >= 0; i--) {
        limb t = sub_mod(x[i], x[n - 1], mods[i]);
        y[i] = mul_mod(t, c[i], mods[i]);
        M[i] = mul_mod(mods[n - 1], c[i], mods[i]);
    }
         
    for(i = n - 3; (int)i >= 0; i--) {
        for(j = i; (int)j >= 0; j--) {
            limb t = mul_mod(y[i + 1], M[j], mods[j]);
            y[j] = sub_mod(y[j], t, mods[j]);
            M[j] = mul_mod(M[j], mods[i + 1], mods[j]);
        }
    }
*/

}


#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template < >
void print_vector< CORE::BigInt >(const std::vector< CORE::BigInt >& v,
        bool hex) {

    unsigned n = v.size();
    printf("(n: %d; ", n);

    for(unsigned j = 0; j < n; j++) 
        printf("%s ", v[j].get_str().c_str());
    printf(")\n");
}
#endif 

template < >
void print_vector< zmod >(const std::vector< zmod >& v, bool no_sign) {

    unsigned n = v.size();

    printf("(n: %d; ", n);
    unsigned m = zmod::modulus;
    for(unsigned j = 0; j < n; j++)
        if(v[j].x < m/2 || no_sign)
            printf("%d ", v[j].x);
        else
            printf("%d ", (v[j].x - m));
    printf(")\n");
}

std::ostream& operator <<(std::ostream& os, const std::vector< zmod >& p) {

    unsigned m = zmod::modulus;
    for(unsigned i = p.size() - 1; i != -1u; i--) {
        zmod t = p[i];
        if(t == zmod(0) && i != 0) 
            continue;
        
        if(i != p.size() - 1) {
            if(t.x < m/2)
                os << " + ";
            else 
                os << " ";
        }
        os << t;

        if(i >= 1) {
            os << "*x";
            if(i > 1)
                os << "^" << i;
        }
    }
    return os;
}

std::ostream& operator <<(std::ostream& os, const zmod& x) {
    
    unsigned m = zmod::modulus;
    if(x.x < m/2)
        os << x.x;
    else
        os << (int)(x.x - m);
    return os;
}

umod_t zmod::modulus(-1u);

// #ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

gmp_randstate_t rands;

// #endif

} // namespace internal

} // namespace CGAL

