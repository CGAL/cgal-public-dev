/*! \file anisotropic_diagram.cpp
 * \brief  anisotropic diagram example.
 */

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Anistropic_voronoi_traits_2.h>

using namespace std;

using AK = CGAL::Arithmetic_kernel;
using Integer = AK::Integer;
using Rational = AK::Rational;

// Definition of Algebraic_kernel_2 (Algebraic_curve_kernel_2)
using Coefficient = Integer;
using Algebraic_curve_kernel_2 =
  CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
  Algebraic_curve_kernel_with_qir_and_bitstream_2;
using Curved_kernel_2 =
  CGAL::Curved_kernel_via_analysis_2<Algebraic_curve_kernel_2>;

using Point_2 = Curved_kernel_2::Point_2;
using Curve_2 = Curved_kernel_2::Curve_2;
using X_monotone_curve_2 = Curved_kernel_2::X_monotone_curve_2;
using VD_traits_2 = CGAL::Anisotropic_voronoi_traits_2<Curved_kernel_2>;
using Anisotropic_diagram_2 = CGAL::Envelope_diagram_2<VD_traits_2>;
using Site_2 = VD_traits::Site_2;

int main(int argc, char* argv[]) {
  std::list<Site_2> sites; // anisotropic sites
  int n;
  std::cin >> n;
  for (int i = 0; i < n; ++i) {
    Raional x, y, a, b, c, d;
    std::cin >> x >> y >> a >> b >> c >> d;
    std::cout << x << " " << y << " " << a << " " << b << " " << c << " " << d
              << endl;

    // making the matrix possitive.
    CGAL_assertion_msg(a*c - b*b > 0, "input invalid");
    sites.push_back(Site_2(make_pair(x, y), a, b, c, d));
  }

  Anisotropic_diagram_2 VD;
  CGAL::lower_envelope_3(sites.begin(), sites.end(), VD);
  std::cout << "Number of sites: " << sites.size() << std::endl;
  std::cout << "Anisotropic VD:" << std::endl
            << "V = " << VD.number_of_vertices()
            << ", E = " << VD.number_of_edges()
            << ", F = " << VD.number_of_faces() << std::endl;

  return 0;
}
