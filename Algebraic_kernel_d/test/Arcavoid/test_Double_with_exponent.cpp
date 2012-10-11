#include <CGAL/basic.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Double_with_exponent.h>

int main () {
  using namespace std;
  using namespace CGAL;

  Double_with_exponent x, y, z;
  x = 2;
  y = 4.1;
  z = 0.0001;

  cout << x << endl;
  cout << y << endl;
  cout << z << endl;
  cout << x*y << endl;
  cout << y.sqrt() << endl;
  cout << z-x << endl;

  for (int i = 0; i < 8; ++i)
    x *= x;
  cout << z+x << endl;

  cout << z/x << endl;
  cout << (-(z-x)).sqrt() << endl;
  x = Double_with_exponent (.5, 258);
  cout << x << endl;
  cout << x.root_d (257) << endl;
  cout << x.root_d (43) << endl;

  leda_integer Iasdf = 1234567;
  Iasdf *= Iasdf * Iasdf;
  Iasdf *= Iasdf * Iasdf;
  leda_bigfloat Fasdf (Iasdf, -1234);
  cout << Fasdf << endl;
  y = Double_with_exponent (Fasdf);
  cout << y << endl;
  Fasdf = y.sqrt();
  cout << Fasdf << endl;

  return 0;
}
