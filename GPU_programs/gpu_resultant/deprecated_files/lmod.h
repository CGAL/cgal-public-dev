
#ifndef LMOD_H
#define LMOD_H

#include <stdio.h>
#include "macros.h"

#define NBITS (sizeof(limb)*8) 

inline void copylimb(limb *dst, const limb *src, int n) {

    memcpy(dst, src, n*sizeof(limb));
}

inline void zerolimb(limb *dst, int n) {

    memset(dst, 0, n*sizeof(limb));
}

// adds two limb vectors modulo 2^(n*NBITS) - 1
// in-place modificaion is allowed
static inline void
add_n_modF(limb *res, const limb *a, const limb *b, int n) {

    limb cy = add_n(res, a, b, n);
    if(cy > 0)
        add_1(res, res, cy, n);
}

// subtracts two limb vectors modulo 2^(n*NBITS) - 1
// in-place modificaion is allowed
static inline void
sub_n_modF(limb *res, const limb *a, const limb *b, int n) {

    limb cy = sub_n(res, a, b, n);
    if(cy > 0)
        sub_1(res, res, cy, n);
}

static inline void
mul_n_modF(limb *res, const limb *a, const limb *b, int n) {

}

//! class storing variable-sized residue modulo \c 2^(32*SZ)-1
//! in semi-normalized form
template < int SZ >
struct lmod { 

    static const int n_words = SZ;

    explicit lmod() { }

    explicit lmod(const limb* _v) {
        copylimb(v, _v, SZ);
    }

    explicit lmod(limb _x) {
        zerolimb(v, SZ);
        v[0] = _x;
    }

    lmod(const lmod& m) {
        copylimb(v, m.v, SZ);
    }

    friend inline lmod & operator += (lmod& z, const lmod& h) { 
        z = z + h; 
        return z;
    }

    friend inline lmod & operator -= (lmod &z, const lmod &h) {
        z = z - h;
        return z;
    }             
 
    friend inline lmod & operator *= (lmod &z, const lmod &h) { 
        z = z * h; 
        return z;    
    }

    friend inline lmod operator - (const lmod &h1, const lmod &h2) { 
        lmod z;
        sub_n_modF(z.v, h1.v, h2.v, SZ); 
        return z;
    }
    
    friend inline bool operator == (const lmod &h1, const lmod &h2) { 

        return (memcmp(&h1.v, &h2.v, SZ*sizeof(limb)) == 0);
    }

    friend inline bool operator != (const lmod &h1, const lmod &h2) { 

        return !(h1 == h2);
    }

    friend inline lmod operator + (const lmod &h1, const lmod &h2) { 
        lmod z;
        add_n_modF(z.v, h1.v, h2.v, SZ); 
        return z;
    }

    friend inline lmod operator * (const lmod &h1, const lmod &h2) { 
        lmod z;
        mul_n_modF(z.v, h1.v, h2.v, SZ); 
        return z;
    }

    lmod & negate() { 
        v = sub_n_modF(lmod(0u), v, SZ);
        return *this; 
    }

    friend inline lmod operator - (const lmod &h) { 
        lmod n(h);  
        n.negate();  
        return n; 
    }

    friend inline lmod operator + (const lmod &h)  { return h; }

    void format(char *buf) {
        char *p = buf;
        for(int i = 0; i < SZ; i++) {
            p += sprintf(p, "0x%x ", v[i]);
        }
        *p = 0;
    }

//private:
    limb v[SZ];   // SZ words to represent residue in semi-normalized form
};

#endif