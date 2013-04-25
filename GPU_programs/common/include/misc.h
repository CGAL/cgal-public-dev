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
    
#ifndef _MISC_H_
#define _MISC_H_

namespace CGAL {

namespace internal {

//! compares 2D arrays of data, \c n elements per row stored with \c padding
//! number of rows given by \c n_batches
//! \c print_when_differs :  indicates whether print elements only if they
//! differ (default)
//! \c print_max : maximal # of entries to print
template < class NT >
bool checkme(const NT *checkit, const NT *truth, unsigned n, unsigned padding,
        unsigned n_batches, NT modulus = 0, unsigned print_max = -1u, 
            bool print_when_differs = true) {

    if(checkit == NULL)
        return false;

    bool res = true;
    unsigned printed = 0;
    printf("\nlegend: batch_id (element_in_batch)\n");
    for(int j = n_batches - 1; j >= 0; j--) {

        const NT *pcheckit = checkit + j * padding,
            *ptruth = truth + j * padding;
        for(int i = n - 1; i >= 0; i--) {

            NT diff = pcheckit[i] - ptruth[i];
            if(diff != 0)
                res = false;
            
// NOTE: %#x - prints 0x[hexval]
            if((diff != 0 || !print_when_differs) && printed < print_max) {
                NT check = pcheckit[i], truth = ptruth[i];
                if(modulus != 0) {
                    if(check >= modulus/2)
                        check -= modulus;
                    if(truth >= modulus/2)
                        truth -= modulus;
                }

                printed++;
                printf("%d (%d) (GPU, truth): %d and %d; diff: %x ", j, i,
                    check, truth, diff);
                if(diff != 0)
                    printf(" DIFFERS\n");
                else 
                    printf("\n");
            }
        }
    }
    return res;
}

template < class NT >
void writeout(const std::vector< NT >& a, unsigned *out, unsigned stride,
         unsigned i_start, unsigned i_beyond, NT div = NT(1)) {

    if(out == 0)
        return;
    if(div != NT(1)) {
        div = NT(mod_inverse(div.x, zmod::modulus));
    }

    for(unsigned i = i_start; (int)i < (int)i_beyond; i++, out += (int)stride) {
        out[0] = (a[i] * div).x;
    }
}

template < class NT >
void fillout(const NT& a, unsigned *out, unsigned stride, unsigned n) {

    if(out == 0)
        return;
    for(unsigned i = 0; (int)i < (int)n; i++, out += stride) {
        out[0] = a.x;
    }
}

#if 0
//! writes \c n (at most 4096) plain words of \c data into file
void save_testcase(const char *filename, const unsigned *data, unsigned n) {

    FILE *fp = fopen(filename, "wb");
    if(fp == NULL) {
        printf("ERROR: unable to open %s..\n", filename);
        return;
    }

    n = (n > 4096 ? 4096 : n);
    fprintf(fp, "%d\n", n);

    for(int i = 0; i < n; i++) {
        fprintf(fp, "%d\n", data[i]);
    }
    fclose(fp);
}

//! reads out plain data from file: provided that  \c data has enough space
//! for \c n elements, returns # of elements being read
void load_testcase(const char *filename, unsigned *data, unsigned n) {

    FILE *fp = fopen(filename, "rb");
    if(fp == NULL) {
        printf("ERROR: unable to open %s..\n", filename);
        return;
    }

    unsigned m;
    fscanf(fp, "%d\n", &m);
    if(m > n)
        m = n;

    for(int i = 0; i < m; i++) {
        fscanf(fp, "%d\n", data + i);
    }
    fclose(fp);
}
#endif

} // namespace internal

} // namespace CGAL

#endif // _MISC_H_
