#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <crlibm.h>

int main (int argc, char **argv) {
  double a, a_10_rndd, a_10_rndn, a_10_rndu;
  double b, b_1_10_rndd, b_1_10_rndn, b_1_10_rndu;

  crlibm_init();

  a = 1.234;
  b = 1.234;

  a_10_rndn = pow_rn (a, 10.);
  b_1_10_rndn = pow_rn (b, .1);

  fesetround (FE_DOWNWARD);

  a_10_rndd = pow (a, 10.);
  b_1_10_rndd = pow (b, .1);

  fesetround (FE_UPWARD);

  a_10_rndu = pow_rn (a, 10.);
  b_1_10_rndu = pow_rn (b, .1);

  printf ("%f ^ 10: %f <= %f ? %i \t <= %f ? %i\n",
          a, a_10_rndd, a_10_rndn, (a_10_rndd <= a_10_rndn),
          a_10_rndu, (a_10_rndn <= a_10_rndu));
  printf ("%f ^ 10: %f <= %f ? %i \t <= %f ? %i\n",
          b, b_1_10_rndd, b_1_10_rndn, (b_1_10_rndd <= b_1_10_rndn),
          b_1_10_rndu, (b_1_10_rndn <= b_1_10_rndu));

  return 0;
}
