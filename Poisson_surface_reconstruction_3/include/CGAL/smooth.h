#include <cmath>
#include <cassert>

#include <array>
namespace CGAL{
/// map multi-index symbols (i,j,k,p), i+j+k+p=3, to linear indices
 enum MultiIndex3 {
   I3000 = 0, I2100 = 1, I1200 = 2, I0300 = 3,
   I2010 = 4, I1110 = 5, I0210 = 6,
   I1020 = 7, I0120 = 8,
   I0030 = 9,

   I2001 = 10, I1101 = 11, I0201 = 12,
   I1011 = 13, I0111 = 14,
   I0021 = 15,

   I1002 = 16, I0102 = 17,
   I0012 = 18,

   I0003 = 19,
 };

/// compute `n!`
 static constexpr int factorial(int n) {
   return n <= 1 ? 1 : (n * factorial(n - 1));
 }

/// compute product `i!*j!*k!*p!`
 static constexpr int prodfac(int i,int j,int k,int p) {
   return factorial(i)*factorial(j)*factorial(k)*factorial(p);
 }

 /// functor for computing powers of an integer
 template <typename FT,int N> struct ipow {
   FT operator()(FT w) { return w*ipow<FT,N-1>()(w); }
 };
 template <typename FT> struct ipow<FT,0> {
   FT operator()(FT w) { return FT(1); }
 };

/** Evaluate cubic trivariate Bernstein polynomial associated with
 * coefficient `(I,J,K,P)` at barycentric coordinate `w`.
 */

 template <typename FT,int I,int J,int K,int P>
 struct bernstein3 {
  FT operator() (const FT* w) {
    return FT(6)/FT(prodfac(I,J,K,P)) *
        ipow<FT,I>()(w[0])*ipow<FT,J>()(w[1]) *
        ipow<FT,K>()(w[2])*ipow<FT,P>()(w[3]);
   }
 };


/** Evaluate trivariate cubic Bernstein polynomial `b` at
    barycentric coordinates `w`.
    \param b coefficients (20 components)
    \param w barycentric coordinates (4 components)
 */
 template <typename FT>
 FT eval_bernstein3(const FT* b,const FT* w) {
   return bernstein3<FT,3,0,0,0>()(w)*b[0] +
         bernstein3<FT,2,1,0,0>()(w)*b[1] +
         bernstein3<FT,1,2,0,0>()(w)*b[2] +
         bernstein3<FT,0,3,0,0>()(w)*b[3] +
         bernstein3<FT,2,0,1,0>()(w)*b[4] +
         bernstein3<FT,1,1,1,0>()(w)*b[5] +
         bernstein3<FT,0,2,1,0>()(w)*b[6] +
         bernstein3<FT,1,0,2,0>()(w)*b[7] +
         bernstein3<FT,0,1,2,0>()(w)*b[8] +
         bernstein3<FT,0,0,3,0>()(w)*b[9] +
         bernstein3<FT,2,0,0,1>()(w)*b[10] +
         bernstein3<FT,1,1,0,1>()(w)*b[11] +
         bernstein3<FT,0,2,0,1>()(w)*b[12] +
         bernstein3<FT,1,0,1,1>()(w)*b[13] +
         bernstein3<FT,0,1,1,1>()(w)*b[14] +
         bernstein3<FT,0,0,2,1>()(w)*b[15] +
         bernstein3<FT,1,0,0,2>()(w)*b[16] +
         bernstein3<FT,0,1,0,2>()(w)*b[17] +
         bernstein3<FT,0,0,1,2>()(w)*b[18] +
         bernstein3<FT,0,0,0,3>()(w)*b[19];
 }

/// compute dot product` <a,b>`
 template <typename FT,size_t N>
 FT dot(const FT* a,const std::array<FT,N>& b) {
   FT sum(0);
   for (size_t i=0;i<N;++i)
     sum+=a[i]*b[i];
   return sum;
 }

/** Compute dot product <a[:,j],b> of matrix column and vector.
    \param a pointer to array that represents `3×N` column-major matrix
    (with leading dimension `N`)
    \param j column index
    \param b vector
*/
 template <typename FT,size_t N>
 FT dot(const FT* a,int j,const std::array<FT,N>& b) {
   FT sum(0);
   for (size_t i=0;i<N;++i)
     sum+=a[i+j*N]*b[i];
   return sum;
 }

/** Compute vector of edges `(i,j)` spanned by `x`.
    \param i index of source vertex
    \param j index of destination vertex
    \param x `3×4` matrix (column-major) with vertex positions in columns
 */
 template <typename FT>
 std::array<FT,3> edge_vector(int i,int j,const FT* x) {
   return std::array<FT,3> { { x[0+j*3]-x[0+i*3],
                                 x[1+j*3]-x[1+i*3],
                                 x[2+j*3]-x[2+i*3] } };
 }




/** Compute control points for cubic element.
    \param[out] b control points/coefficients of trivariate cubic
    polynomial (20 components)
    \param x `3×4` column-major matrix spanning tetrahedron (vertex
    positions in columns)
    \param f vector of 4 samples of function values at vertices
    \param gradf `3×4` column-major matrix with samples of gradients
    at vertices in columns
 */
 template <typename FT>
 void control_points(FT* b,const FT* x,const FT* f,const FT* gradf) {
  // End point interpolation.
   b[I3000] = f[0];
   b[I0300] = f[1];
   b[I0030] = f[2];
   b[I0003] = f[3];

  // Interpolate gradient at end points.

   auto e01 = edge_vector(0,1,x);  // edge vectors:
   auto e12 = edge_vector(1,2,x);  //  eij refers to edge i->j
   auto e20 = edge_vector(2,0,x);
   auto e03 = edge_vector(0,3,x);
   auto e13 = edge_vector(1,3,x);
   auto e23 = edge_vector(2,3,x);

   auto d01 = +dot(gradf,0,e01);   // scaled (3*) directional derivatives:
   auto d02 = -dot(gradf,0,e20);   //  dij refers to derivative in direction
   auto d03 = +dot(gradf,0,e03);   //  of edge i->j (eij) evaluated at vertegradf i

   auto d10 = -dot(gradf,1,e01);
   auto d12 = +dot(gradf,1,e12);
   auto d13 = +dot(gradf,1,e13);

   auto d20 = +dot(gradf,2,e20);
   auto d21 = -dot(gradf,2,e12);
   auto d23 = +dot(gradf,2,e23);

   auto d30 = -dot(gradf,3,e03);
   auto d31 = -dot(gradf,3,e13);
   auto d32 = -dot(gradf,3,e23);

   b[I2100] = d01/FT(3)+b[I3000];
   b[I2010] = d02/FT(3)+b[I3000];
   b[I2001] = d03/FT(3)+b[I3000];

   b[I1200] = d10/FT(3)+b[I0300];
   b[I0210] = d12/FT(3)+b[I0300];
   b[I0201] = d13/FT(3)+b[I0300];

   b[I1020] = d20/FT(3)+b[I0030];
   b[I0120] = d21/FT(3)+b[I0030];
   b[I0021] = d23/FT(3)+b[I0030];

   b[I1002] = d30/FT(3)+b[I0003];
   b[I0102] = d31/FT(3)+b[I0003];
   b[I0012] = d32/FT(3)+b[I0003];

   // Determine remaining coefficients such that quadratic polynomials
   // are reproduced.

   auto c1100 = (b[I3000]+(d01+d10)/FT(2)+b[I0300])/FT(2);
   auto c0110 = (b[I0300]+(d12+d21)/FT(2)+b[I0030])/FT(2);
   auto c1010 = (b[I3000]+(d20+d02)/FT(2)+b[I0030])/FT(2);

   auto c1001 = (b[I3000]+(d03+d30)/FT(2)+b[I0003])/FT(2);
   auto c0101 = (b[I0300]+(d13+d31)/FT(2)+b[I0003])/FT(2);
   auto c0011 = (b[I0030]+(d23+d32)/FT(2)+b[I0003])/FT(2);

   b[I1110] = (c1010+c1100+c0110)/FT(3);
   b[I1101] = (c1001+c1100+c0101)/FT(3);
   b[I1011] = (c1001+c1010+c0011)/FT(3);
   b[I0111] = (c0101+c0110+c0011)/FT(3);
 }

/** Evaluate cubic element defined by `,f,gradf` at barycentric coordinates `w`.

    The element interpolates function values `f` and gradients `gradf` at the
    vertices of the tetrahedron spanned by `x`.

    The interpolant is C^1 at vertices and reproduces quadratic polynomials.

    \param x `3×4` column-major matrix spanning tetrahedron (vertex
    positions in columns)
    \param f vector of 4 samples of function values at vertices
    \param gradf `3×4` column-major matrix with samples of gradients
    at vertices in columns
    \param w 4-vector of barycentric coordinates

    \return interpolated function value
 */
 template <typename FT>
 FT eval_element3(const FT* x,const FT* f,const FT* gradf,
                   const FT* w) {
   FT b[20];
   control_points(b,x,f,gradf);
   return eval_bernstein3(b,w);

 }
}
