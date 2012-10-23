

// Vandermonde interpolation inner loop (4 points per thread)
//! set x0 = nn - j !!!!
#if 1
    // problem: large warp divergence with this pattern..
    // in fact only a single thread diverges at a time
    if((int)a0 >= 0) {
        t1 = b0 + 1 - x0;  // 2-way BCs when reading xs's, bad..
        //NOTE NOTE you may also take into account that t1 is typically small
        // i.e., 1 <= t1 <= 10..
        a.x = mul_small_m(a.x, t1, m, invk, mx2, _2e23);
    }

    if((int)a0 >= 1) {
        t1 = b0 + 2 - x0; // we know that evaluation points are
                                // monotonically increasing 
        a.y = mul_small_m(a.y, t1, m, invk, mx2, _2e23);
    }

    if((int)a0 >= 2) {
        t1 = b0 + 3 - x0; // we know that evaluation points are
                                // monotonically increasing 
        a.z = mul_small_m(a.z, t1, m, invk, mx2, _2e23);
    }

    if((int)a0 >= 3) {
        t1 = b0 + 4 - x0; // we know that evaluation points are
                                // monotonically increasing 
        a.w = mul_small_m(a.w, t1, m, invk, mx2, _2e23);
    }

    if(thid == (j >> 2)) {
        a_prev = m;
    }

    if((int)a0 < 0) {
        t1 = a.x; // save a.x
        a.x = mul_small_m(a.x, x0, m, invk, mx2, _2e23);
        a.x = sub_m(a_prev, a.x, m);
        a_prev = t1;
    }

    if((int)a0 <= 0) {
        // a_prev not used anymore
        t1 = a.y;
        a.y = mul_small_m(a.y, x0, m, invk, mx2, _2e23);
        a.y = sub_m(a_prev, a.y, m);
        a_prev = t1;
    }

    if((int)a0 <= 1) {
        t1 = a.z; // t1 not used anymore
        a.z = mul_small_m(a.z, x0, m, invk, mx2, _2e23);
        a.z = sub_m(a_prev, a.z, m);
        a_prev = t1;
    }

    if((int)a0 <= 2) {
        a.w = mul_small_m(a.w, x0, m, invk, mx2, _2e23);
        a.w = sub_m(a_prev, a.w, m);
    }
#endif

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// 512-point interpolation using 256 threads (inner loop):

    if((int)a0 > 0) {
        t1 = b0 + 1 - x0;  // 2-way BCs when reading xs's, bad..
        //NOTE NOTE you may also take into account that t1 is typically small
        // i.e., 1 <= t1 <= 10..
        a.x = mul_small_m(a.x, t1, m, invk, mx2, _2e23);

        t1 = b0 + 2 - x0; // we know that evaluation points are
                                // monotonically increasing 
        a.y = mul_small_m(a.y, t1, m, invk, mx2, _2e23);

    } else if(a0 == 0) {

        t1 = b0 + 1 - x0; // 2-way BCs.. bad ((
        a.x = mul_small_m(a.x, t1, m, invk, mx2, _2e23);

        t1 = a_prev;
        if(thid == (j >> 1)) {
            t1 = m;
        }
        a.y = mul_small_m(a.y, x0, m, invk, mx2, _2e23);
        a.y = sub_m(t1, a.y, m);
    } else {

        if(thid == (j >> 1)) {
            a_prev = m;
        }
        // does the order of evaluation makes sense ??
        // I guess not..
        t1 = a.x; // save a.x
        a.x = mul_small_m(a.x, x0, m, invk, mx2, _2e23);
        a.x = sub_m(a_prev, a.x, m);

        a.y = mul_small_m(a.y, x0, m, invk, mx2, _2e23);
        a.y = sub_m(t1, a.y, m);
    }

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //! this just defines the order of evaluation points - can be changed later
    x0 = j + 1; //nn - j; //xs[j];
    // from here on, xs[j] is not used anymore..

    if(j == thid)
        w.x = a0;
    else if(j == thid + BlockSz)
        w.y = a0;
    else if(j == thid + BlockSz*2)
        w.z = a0;
    else if(j == thid + BlockSz*3)
        w.w = a0;

    // b0 not needed anymore
    limb b0 = __umul24(thid, 0x1000004) + 1;
    a0 = j - b0; // == n_m_j - (thid * 4 + 1)
    b0 += nn - j - 1;

    if(thid == (j >> 2)) {
        a_prev = m;
    }

    j--; // in fact on the last iteration we can skip all these steps
    if(j == -1) //! number of iterations now is WS*16 !!
        break;

/**
 if((int)a0 > 2) { 4 first

    } else if((int)a0 == 2) { // 3 first, 1 second

    } else if((int)a0 == 1) { // 2 first, 2 second

    } else if((int)a0 == 0) { // 1 first, 3 second

    a0 < 0: 4 second
*/

    limb t1 = m, t2 = m, t3 = m, t4 = m,
        mul1 = x0, mul2 = x0, mul3 = x0, mul4 = x0;

    if((int)a0 >= 0) {

        // TODO: any way to replace this by __sad ??
        // mul[i] is either x0 or x0 - x[i]
        mul1 -= (N - b0);
    } //TODO: use else here..

    // we do not need b.x anymore => used as a new a.x
    Bs[0] = b.x; 
    b.x = mul_small_m(mul1, a.x, m, invk, mx2);

    mul1 = m; // use mul1 instead of t1

    if((int)a0 >= 1) {
        mul2 -= (N - 1 - b0);
    }

    if((int)a0 >= 2) {
        mul3 -= (N - 2 - b0);
    }

    if((int)a0 >= 3) {
        mul4 -= (N - 3 - b0);
    }
    // b0 is a new a.y here..
    b0 = mul_small_m(mul2, a.y, m, invk, mx2);

    mul2 = m; // reuse mul2
    if((int)a0 < 0) {
        mul1 = a_prev;
        mul2 = a.x;
    }

    if((int)a0 == 0) {
        mul2 = a_prev;
    }

    if((int)a0 < 1) {
        t3 = a.y;
    }    

    a.x = sub_m(mul1, b.x, m); // b.x is a new a.x

    // at this point we can per

    if((int)a0 == 2) {
        t4 = a_prev;
    }

    if((int)a0 < 2) {
        t4 = a.z;
    }

    a.y = sub_m(mul2, b0, m);

    if((int)a0 == 1) {
        t3 = a_prev;
    }

    a.z = mul_small_m(mul3, a.z, m, invk, mx2);
    a.z = sub_m(t3, a.z, m);

    a.w = mul_small_m(mul4, a.w, m, invk, mx2);
    a.w = sub_m(t4, a.w, m);

    //! imagine you have to perform all these instructions at *each* iteration
    //! does it cost one mul_m ?? (taking into account that you will have
    //! to execute reduction with 256 threads at the end..

    // NOTE: possible read-write hazard on the next iteration
    As[0] = a.x, a_prev = a.x;

    CU_SYNC

    // shift down the variables
    a.x = a.y, b.x = b.y, a.y = a.z, b.y = b.z;
    a.z = a.w, b.z = b.w;

        //TODO: can use a0 & b0 here for temporary variables..
    if(thid + 1 < n_thids) { // TODO: can merge these two conditions into one
        a.w = As[1]; // read out elements shifted by one (except the last
        b.w = Bs[1]; // thid)
    }

    a0 = r[BlockSz]; // == As[-thid]

    } // while(1)

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inner loop original:

 while(1) {

    if(thid == n_thids - 1) {
        a.w = 1;
        b.w = m; // so we have 0 * a0 - 1 * b0 = -b0
    }

    //! shift all index ranges by -(j+1)
    // r[0] == b0
    limb x0 = r[0];

    limb mx100 = __umul24(100, m), mx2 = 0;//__umul24(0x1000002, m);
    b.x = sub_mul_reduce_m(b.x, a0, a.x, x0, m, invk, invm, mx100);
    b.y = sub_mul_reduce_m(b.y, a0, a.y, x0, m, invk, invm, mx100);
    b.z = sub_mul_reduce_m(b.z, a0, a.z, x0, m, invk, invm, mx100);
    b.w = sub_mul_reduce_m(b.w, a0, a.w, x0, m, invk, invm, mx100);

    //! this just defines the order of evaluation points - can be changed later
    x0 = j + 1; //nn - j; //xs[j];
    // from here on, xs[j] is not used anymore..

    if(j == thid)
        w.x = a0;
    else if(j == thid + BlockSz)
        w.y = a0;
    else if(j == thid + BlockSz*2)
        w.z = a0;
    else if(j == thid + BlockSz*3)
        w.w = a0;

    // b0 not needed anymore
    limb b0 = __umul24(thid, 0x1000004) + 1;
    a0 = j - b0; // == n_m_j - (thid * 4 + 1)
    b0 += nn - j - 1; 

    if(thid == (j >> 2)) {
        a_prev = m;
    }

/**
 if((int)a0 > 2) { 4 first

    } else if((int)a0 == 2) { // 3 first, 1 second

    } else if((int)a0 == 1) { // 2 first, 2 second

    } else if((int)a0 == 0) { // 1 first, 3 second

    a0 < 0: 4 second
*/

    limb t1 = m, t2 = m, t3 = m, t4 = m,
        mul1 = x0, mul2 = x0, mul3 = x0, mul4 = x0;

    if((int)a0 >= 0) {

        // TODO: any way to replace this by __sad ??
        // mul[i] is either x0 or x0 - x[i]
        mul1 -= (N - b0);
    } //TODO: use else here..

    if((int)a0 >= 1) {
        mul2 -= (N - 1 - b0);
    }

    if((int)a0 >= 2) {
        mul3 -= (N - 2 - b0);
    }

    if((int)a0 >= 3) {
        mul4 -= (N - 3 - b0);
    }

    // at this point we can per

    if((int)a0 == 2) {
        t4 = a_prev;
    }

    if((int)a0 < 2) {
        t4 = a.z;
    }

    if((int)a0 < 1) {
        t3 = a.y;
    }    

    if((int)a0 == 1) {
        t3 = a_prev;
    }

    if((int)a0 == 0) { 
        t2 = a_prev;
    }

    if((int)a0 < 0) {
        t1 = a_prev;
        t2 = a.x;
    }

    a.x = mul_small_m(mul1, a.x, m, invk, mx2);
    a.x = sub_m(t1, a.x, m);

    a.y = mul_small_m(mul2, a.y, m, invk, mx2);
    a.y = sub_m(t2, a.y, m);

    a.z = mul_small_m(mul3, a.z, m, invk, mx2);
    a.z = sub_m(t3, a.z, m);

    a.w = mul_small_m(mul4, a.w, m, invk, mx2);
    a.w = sub_m(t4, a.w, m);

    //! imagine you have to perform all these instructions at *each* iteration
    //! does it cost one mul_m ?? (taking into account that you will have
    //! to execute reduction with 256 threads at the end..

    j--;
    if(j == -1) //! number of iterations now is WS*16 !!
        break;
    
    // NOTE: possible read-write hazard on the next iteration
    As[0] = a.x, Bs[0] = b.x, a_prev = a.x;

    CU_SYNC

    // shift down the variables
    a.x = a.y, b.x = b.y, a.y = a.z, b.y = b.z;
    a.z = a.w, b.z = b.w;

        //TODO: can use a0 & b0 here for temporary variables..
    if(thid + 1 < n_thids) { // TODO: can merge these two conditions into one
        a.w = As[1]; // read out elements shifted by one (except the last
        b.w = Bs[1]; // thid)
    }

    a0 = r[BlockSz]; // == As[-thid]

    } // while(1)



 template< class Coeff >
    struct Simple_matrix
        : public std::vector< std::vector< Coeff > > {
        
        typedef Coeff NT;    

        Simple_matrix() {
        }
                
        Simple_matrix( int m ) {
            initialize( m, m, Coeff(0) );
        }
                        
        Simple_matrix( int m, int n, Coeff x = Coeff(0) ) {
            initialize( m, n, x );
        }
        
        void swap_rows(int i, int j) {
            std::vector< Coeff > swap = this->operator[](i);
            this->operator[](i) = this->operator[](j);
            this->operator[](j) = swap;
        }

        void swap_columns(int i, int j) {
            for(int k = 0; k < m; k++) {
                Coeff swap = this->operator[](k).operator[](i);
                this->operator[](k).operator[](i)
                    = this->operator[](k).operator[](j);
                this->operator[](k).operator[](j) = swap;
            }
        }

        //! row - outer dimension
        //! column - inner dimension
        //!: vector< vector < >(n_cols) > (n_rows)

        unsigned n_rows() const { return m; }
        unsigned n_cols() const { return n; }
        
        private:
        
        void initialize( int m, int n, Coeff x ) {
            this->reserve( m );
            this->m = m;
            this->n = n;
            for( int i = 0; i < m; ++i ) {
                this->push_back( std::vector< Coeff >() );
                this->operator[](i).reserve(n);
                for( int j = 0; j < n; ++j ) {
                    this->operator[](i).push_back( x );
                }
            }            
        }
        
        int m, n;
    };

    template< class Coeff >
    struct Simple_vector
        : public std::vector< Coeff > {

        typedef Coeff NT;

        Simple_vector() {
        }

        Simple_vector(int m) : 
             std::vector< Coeff >(m) {
        }

        Simple_vector(int m, Coeff c) : 
             std::vector< Coeff >(m, c) {
        }

        int degree() const {
//             if(this->size() == 0)
//                 return 0;
            return this->size() - 1;
        }

        Coeff operator*( const Simple_vector<Coeff>& v2 ) const {
            Coeff result(0);
                        
            for( unsigned i = 0; i < this->size(); ++i )
                result += ( this->operator[](i) * v2[i] );
                
            return result;
        }
    };

template < class NT >
void print_matrix(const Simple_matrix< NT >& mt) {

    unsigned m = mt.n_rows(), n = mt.n_cols();
    printf("============ nrows: %d; ncols: %d ============\n", m, n);

    for(unsigned i = 0; i < m; i++) {

        for(unsigned j = 0; j < n; j++) {
            printf("%d ", mt[i][j]);
        }
        printf("\n");
    }
    printf("========================\n");
}


template < >
void print_matrix(const Simple_matrix< Integer >& mt) {

    unsigned m = mt.n_rows(), n = mt.n_cols();
    printf("============ nrows: %d; ncols: %d ============\n", m, n);

    for(unsigned i = 0; i < m; i++) {

        for(unsigned j = 0; j < n; j++) {
            printf("%s ", mt[i][j].get_str().c_str());
        }
        printf("\n");
    }
    printf("========================\n");
}

// Part of det_berkowitz
// Computes sum of all clows of length k
template <class M>
std::vector<typename M::NT> clow_lengths (const M& A, int k, int n) {

    typedef typename M::NT NT;
    
    int i, j, l;
    
    Simple_vector<NT> r(k-1), s(k-1), t(k-1);
    std::vector<NT> rMks(k);
    
    Simple_matrix<NT> MM(k-1);
    
    for(i= n - k + 2; i <= n; ++i)
        for(j = n - k + 2; j <= n; ++j)
            MM[i-n+k-2][j-n+k-2] = A[i-1][j-1];
    
    i = n - k + 1;
    l = 1;
    
    for(j = n - k + 2; j <= n; ++j, ++l) {
        r[l-1] = A[i-1][j-1];
        s[l-1] = A[j-1][i-1];
    }
    
    rMks[0] = A[i-1][i-1];
    rMks[1] = r*s;
    
    for(i = 2; i < k; ++i) {
        // r = r * M;
        for(j=0;j<k-1;++j)
            for(l=0;l<k-1;++l)
                t[j] += r[l] * MM[l][j];

        for(j = 0; j < k - 1; ++j) {
            r[j] = t[j];
            t[j] = NT(0);
        }
        rMks[i] = r*s;
    }
    return rMks;
}

template <class M, class OutputIterator>
OutputIterator minors_berkowitz (const M& A, OutputIterator minors,
        int n, int m = 0)  {

    typedef typename M::NT NT;
        
    // If default value is set, reset it to the second parameter
    if(m == 0) 
        m = n;
        
    int i, j, k, offset;
    std::vector<NT> rMks;
    NT a;
    
    Simple_matrix<NT> B(n+1);
    Simple_vector<NT> p(n+1), q(n+1);
    
    for(k = 1; k <= n; ++k) {
        // compute vector q = B*p;
        if (k == 1) {
                p[0] = NT(0) - NT(1);
                q[0] = p[0];
                p[1] = A[n-1][n-1];
                q[1] = p[1];
        } else if (k == 2) {
                p[0] = NT(1);
                q[0] = p[0];
                p[1] = -A[n-2][n-2] - A[n-1][n-1];
                q[1] = p[1];
                p[2] = -A[n-2][n-1] * A[n-1][n-2] + A[n-2][n-2] * A[n-1][n-1];
                q[2] = p[2];
        } else if (k == n) {
            rMks = clow_lengths<M>(A,k,n);
                // Setup for last row of matrix B
                i = n+1;
                B[i-1][n-1] = NT(0) - NT(1);
    
                for (j=1;j<=n;++j)
                    B[i-1][i-j-1] = rMks[j-1];
    
                p[i-1] = NT(0);
    
                for (j=1;j<=n;++j)
                    p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];
        }
            else
            {
                rMks = clow_lengths<M>(A,k,n);
    
                // Setup for matrix B (diagonal after diagonal)
    
                for (i=1;i<=k;++i)
                    B[i-1][i-1] = NT(0) - NT(1);
    
                for (offset=1;offset<=k;++offset)
                {
                    a = rMks[offset-1]; 
    
                    for (i=1;i<=k-offset+1;++i)
                        B[offset+i-1][i-1] = a;
                }
    
                // Multiply s.t. p=B*q
    
                for (i=1;i<=k;++i)
                {
                    p[i-1] = NT(0);
    
                    for (j=1;j<=i;++j)
                        p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];
                }
    
                p[i-1] = NT(0);
    
                for (j=1;j<=k;++j)
                    p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];
    
                for (i=1;i<=k+1;++i)
                    q[i-1] = p[i-1];
            }
           
        if(k > n-m) { 
          (*minors)=p[k];
          ++minors;
        }
        }
    return minors;
}
    
template <class M> inline
typename M::NT det_berkowitz (const M& A) {

    typedef typename M::NT NT;
    if(A.n_cols() == 0) 
        return NT(1);
      
    NT det[1];
    minors_berkowitz(A, det, A.n_cols(), 1);
    return det[0];
}

template < class NT >
Simple_matrix< NT > hybrid_bezout_matrix(Simple_vector< NT > f,
         Simple_vector< NT > g, int sub = 0) {

    typedef Simple_matrix<NT> Matrix;

    int n = f.size() - 1, m = g.size() - 1;
    int i, j, k, l;

    if(m > n) {
        std::swap(f, g);
        std::swap(m, n);
    }

    Matrix B(n-sub);
    for (i = 1+sub; i <= m; i++) {
        for (j = 1; j <= n-sub; j++) {
            NT s(0);
            for (k = 0; k <= i-1; k++) {
                l = n+i-j-k;
                if ((l <= n) and (l >= n-(m-i))) {
                    s += f[l] * g[k];
                }
            }
            for (k = 0; k <= n-(m-i+1); k++) {
                l = n+i-j-k;
                if ((l <= m) and (l >= i)) {
                    s -= f[k] * g[l];
                }
            }
            B[i-sub-1][j-1] = s;
        }
    }

    for (i = std::max(m+1, 1+sub); i <= n; i++) {
        for (j = i-m; j <= std::min(i, n-sub); j++) {
            B[i-sub-1][j-1] = g[i-j];
        }
    }

    return B; // g++ 3.1+ has NRVO, so this should not be too expensive
}


#if CUMP_COMPILE_TEST_KERNEL

#define UINT2FLOAT(x,mode) \
    XXX::__internal_uint2float_kernel((x), (cudaRoundMode)(mode))

#define FLOAT2UINT(x,mode) \
    XXX::__internal_float2uint((x), (cudaRoundMode)(mode))

#define MUL_ROUND(a, b, rndNearest) \
    XXX::__internal_fmul_kernel(a, b, rndNearest)

unsigned tst_mul(unsigned a, unsigned b, unsigned m,
        unsigned i2f, unsigned f2i)
{
    if(a >= m || b >= m) {
        printf("mul_mod error: out of bounds: a: %u; b: %u; m: %u\n", a, b, m);
        exit(1);
    }

    unsigned mullo = my_umul24(a, b), mulhi = my_umul24hi(a, b), l;

    unsigned k = (m >> 9);
    float invk = 128.0f / k, safe = 1;
    float prodf = UINT2FLOAT(mulhi, i2f);

    prodf = MUL_ROUND(prodf, invk, 0);
    l = FLOAT2UINT(prodf - 0.5f, f2i);

    unsigned truth = (((umod_t )mulhi) << 7) / (umod_t)k;
//  diff = (((umod_t )mulhi) << 7) - (umod_t)k * (umod_t)l;
    unsigned diff = truth - l;
    
  if(truth != l) {
        printf("prodf = %d --- ",(unsigned&)safe);
   printf("a: %x b: %x; prod=%llx truth = %x cdiv2k = %x diff= %d m = %x k: %x\n",
        a, b,  (umod_t)a*(umod_t)b,
            truth, l, diff, m, k);
   }

   return diff;

}

const char *_mode2str(unsigned md) {
    switch(md) {
    case cudaRoundNearest:
        return "rn";
    case cudaRoundZero:
        return "rz";
    case cudaRoundPosInf:
        return "ru";
    case cudaRoundMinInf:
        return "rd";
    }
    return 0;
}

void check_mul_rounding(unsigned a, unsigned b, unsigned m) {

    printf("computing %d * %d mod %x\n", a, b, m);

    for(unsigned i2f = cudaRoundNearest; i2f <= cudaRoundMinInf; i2f++) {
        for(unsigned f2i = cudaRoundNearest; f2i <= cudaRoundMinInf; f2i++) {
            
            printf("\ni2f: %s; f2i %s; res: %d ============= ",
                _mode2str(i2f), _mode2str(f2i),
                    tst_mul(a, b, m, i2f, f2i));
            printf("\n\n==========================");
        }
    }

}
#endif


unsigned mod_inverse_montgomery(unsigned x, unsigned m, unsigned& niters) {

    int u = m, v = x, r = 0, s = 1, k = 0;

    while(v > 0) {

        if((u & 1) == 0) {
            u >>= 1;
            s <<= 1;
        } else if((v & 1) == 0) {
            v >>= 1;
            r <<= 1;
        } else {
            // in case dif == 1 we have u = 0
            int dif = (u - v) >> 1, sum = r + s;
            if(dif > 0) {
                u = dif, r = sum;
                s <<= 1;
            } else {
                v = -dif, s = sum;
                r <<= 1;
            }
        }
        k++;
    }

    niters = k;

    if(r >= m)
        r = r - m;
    r = m - r;

    return r;

#if 0
    unsigned beta = 1 << 24;

    unsigned u1, u2;
    egcd(beta, m, u1, u2);
    if((int)u2 < 0)
        u2 += beta;

    unsigned mu = beta - u2;
    unsigned test = ((uint64)mu * (uint64)m) & 0xffffff;

    if(k > 24) {
        //r = montgomery_mul(r, 1, m, mu);

        unsigned C = r;
        unsigned Q = my_umul24(C, mu);

        limb hi, lo;
        mul_24_wide_host(hi, lo, Q, m);
        lo = my_umul24(lo, 1) + C;
        hi = (hi >> 8) + (lo >= beta);

        //hi >>= 16, hi += (lo < C);
        //hi = (hi << 8) + (lo >> 24); // remap to 24-24
        r = hi;
        k -= 24;
    }
    //r = montgomery_mul(r, 1 << (24 - k), m, mu);

    unsigned c_small = r << (24 - k);
    unsigned Q = my_umul24(c_small, mu);

    limb hi1, lo1, hi2, lo2;
    mul_24_wide_host(hi1, lo1, Q, m);

    lo2 = c_small, hi2 = r >> k; // == c_small ??

    lo1 &= 0xffffff;
    unsigned sum = lo1 + (lo2 & 0xffffff);
    hi1 = (hi1 >> 8) + hi2 + (sum >> 24);

    r = hi1;

    if(r < 0 || r >= m) {
        printf("wrong r: %d\n", r); exit(1); }

    // can iterate until k - 24 then, multiply by precomputed inverse
//     for(int i = 0; i < (int)(k); i++) {
//         if((r & 1) == 0)
//             r >>= 1;
//         else
//             r = (r + m) >> 1;
//     }
    return r;
#endif
}


unsigned mod_inverse_lshift(unsigned x, unsigned m, unsigned& n_iters) {

    unsigned u = 0, v = 0;
    int U = m, V = x, r = 0, s = 1;
    unsigned n = 24;

    n_iters = 0;
    while(std::abs((int)U) != (1 << u) && std::abs((int)V) != (1 << v)) {

        if(std::abs((int)U) < (1 << (n-1))) { // 0x80000000
            U <<= 1, u++;
            if(u > v)
                r <<= 1;
            else
                s >>= 1; //! shound this be a signed value ??

        } else if(std::abs((int)V) < (1 << (n-1))) {
            V <<= 1, v++;
            if(v > u)
                s <<= 1;
            else
                r >>= 1;
        } else {

            if((U & 0x80000000) == (V & 0x80000000)) {

                unsigned duv = U - V, drs = r - s;
                if(u <= v) {
                    U = duv, r = drs;
                } else {
                    V = -duv, s = -drs;
                }
            } else {
                unsigned duv = U + V, drs = r + s;
                if(u <= v) {
                    U = duv, r = drs;
                } else {
                    V = duv, s = drs;
                }
            }
        }
        if(U == 0 || V == 0) return 0;
        n_iters++;
//         printf("U = %x V = %x r = %d s = %d\n", U, V, r, s);
    }

    if(std::abs((int)V) == (1 << v)) {
        r = s, U = V;
    }
    if((int)U < 0) {
        if(r < 0)
            r = -r;
        else
            r = m - r;
    }
    if((int)r < 0)
        r = r + m;
    return (unsigned)r;
}

unsigned mod_inverse_rshift(unsigned x, unsigned m, unsigned& niters) {

    unsigned u = m, v = x, r = 0, s = 1; // TODO: change to int ??
    niters = 0;
    while(v > 0) {

        if((u & 1) == 0) {
            u >>= 1;
            if(r & 1)
                r = (r + m) >> 1;
            else
                r >>= 1;
        } else if((v & 1) == 0) {
            v >>= 1;
            if(s & 1)
                s = (s + m) >> 1;
            else
                s >>= 1;
        } else {
            int dif = u - v;
            if(dif > 0) {
                u = dif, r = r - s;
                if((int)r < 0)
                    r = r + m;
            } else {
                v = -dif, s = s - r;
                if((int)s < 0)
                    s = s + m;
            }
        }
        niters++;
    }

    if(r >= m)
        r = r - m;
    if((int)r < 0)
        r = r + m;
    return r;
}

// mu = -N^-1 mod beta
unsigned montgomery_mul(unsigned x, unsigned y, unsigned N, unsigned mu) {

    uint64 C = (uint64)x * (uint64)y;
    uint64 Q = (C * (uint64)mu) & 0xffffff;
    uint64 R = (Q * (uint64)N + C) >> 24;

    if(R >= (1 << 24))
        R -= N;
    return (unsigned)R;
}

// Set q,t so that  n == q * 2^t + 1
// n must not equal 1, else routine loops.
void n2qt(const umod_t n, umod_t &q, int &t) {
    q = n - 1;  t = 0;
    while ( 0==(q & 1) )  { q >>= 1; ++t; }
}

// Return whether n is a strong pseudoprime to base a.
// q and t must be set so that  n == q * 2^t + 1
bool is_strong_pseudo_prime(const umod_t n, const umod_t a, const umod_t q, 
        const int t) {
    umod_t b = pow_mod(a, q, n);

    if ( 1==b )  return true;  // passed
    // if ( n-1==b )  return true;  // passed

    int e = 1;
    while ( (b!=1) && (b!=(n-1)) && (e<t) )
    {
        b = mul_mod(b, b, n);
        e++;
    }
    if ( b!=(n-1) )  return false;  // =--> composite
    return  true;  // passed
}

// Rabin-Miller compositeness test.
// Return true of none of the bases <=cm prove compositeness.
// If false is returned then n is proven composite (also for n=1 or n=0).
// If true is returned the probability
//   that n is composite is less than (1/4)^cm
bool rabin_miller(umod_t n, uint cm) {

    if ( n<=1 )  return false;

    if((n & 1) == 0)
        return (n == 2);

    umod_t nh = n / 2;
    if(nh < (umod_t)oddprime_array.size())
        return oddprime_array[(ulong)nh];

    //if ( n < small_prime_limit )  return  is_small_prime( (ulong)n );

    umod_t q;
    int t;
    n2qt(n, q, t);

    if ( 0==cm )  cm = 20;  // default
    uint c = 0;
    while ( ++c<=cm )
    {
        umod_t a = c + 1;

        // if n is a c-SPP then it also is a c**k (k>1) SPP.
        // That is, powers of a non-witness are non-witnesses.
        // So we skip perfect powers:
        //if ( is_small_perfpow(a) )  continue;

        if ( a >= n )  return  true;
        if ( !is_strong_pseudo_prime(n, a, q, t) )  return false;  
    }

    return true;  // strong pseudoprime for all tested bases
}

//! /////////////////////////////////////////////////////////////////////////////
//! part of resultant debugging code
//! /////////////////////////////////////////////////////////////////////////////
        mpz_t t;
        mpz_init(t);

        Vector_1 v1(v1m.size()), v2(v2m.size());
        unsigned k;

        for(k = 0; k < v1.size(); k++)
            v1[k] = Integer(v1m[k].x);

        for(k = 0; k < v2.size(); k++)
            v2[k] = Integer(v2m[k].x);

        Matrix bz = hybrid_bezout_matrix(v1, v2);
        Integer truth = det_berkowitz(bz);

        mpz_mod_ui(t, truth.get_mp(), m);
        zmod truth_m(t->_mp_d[0]);

        //if(det != truth_m)
            std::cout << "wrong determinant: " << det << "; truth: " <<
                truth_m << "\n";
            mpz_clear(t);
//! /////////////////////////////////////////////////////////////////////////////
//! part of resultant debugging code
//! /////////////////////////////////////////////////////////////////////////////
