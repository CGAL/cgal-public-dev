

#include <stdio.h>
#include <math.h>

#define DO_NOT_USE_THIS_STUFF
#include "include/ieee754_fp.h"

void print_bin(unsigned x) {

    if(x == 0) 
        printf("0");

    printf("hex = %x: ; bin: ", x);

    while(x) {
        if(x & (1<<31))
            printf("1");
        else
            printf("0");
        x <<= 1;
    }
    printf("\n");
}

void float2int() {

    float _1e23 = (float)(1<<23);
    float lim = 1.0f * _1e23;
    float x =  7777.7676;//_1e23 - 222;

    float y = x + lim;
    int xi = (int&)y;

    printf("lim = "); print_bin((unsigned&)lim);
    printf("x = "); print_bin((unsigned&)x);
    printf("y = "); print_bin((unsigned&)y);

    //xi &= 0x3FFFFF; // mask 23 bits, not 24 !!

    xi = xi - (int&)lim;

//     if(x >= lim) {
//         printf("!!!!!!!!!!!!!overflow\n");
//         xi = (((int&)x) & 0xFFFFFF) + (1<<23);
//     }

    printf("xi = "); print_bin((unsigned&)xi);

    printf("lim = %f x = %f xi = %d \n", lim, x, xi);
    printf("diff = %d\n", (int)floor(x) - xi);
}

void int2float() {
    unsigned x = (1 << 31) + 7771263, xi = x;

    float xf = (float)x;

    //shift = __internal_normalize((unsigned int*)&res.i);

//   t = res.i << 24;
//   res.i = (res.i >> 8);
//   res.i += (127 + 30 - shift) << 23;
    xi = (xi >> 8);
    xi += (127 + 30) << 23;
    float xif = (float&)xi;

    printf("xf = %f x = %u xif = %f\n", xf, x, xif);
}

// TODO: fast float2int & int2float ?
int main() {

    modifyFPUStateX86(__FPU_CW_ROUND_MASK__, __FPU_CW_ROUND_CHOP__);
    int2float();

    return 1;
}
