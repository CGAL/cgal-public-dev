#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
//typedef CGAL::Polyhedron_3<Kernel> Mesh;

typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
typedef boost::property_map<Mesh, Vertex_distance_tag >::type Vertex_distance_map;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Heat_method_3::Heat_method_3<Mesh,Kernel,Vertex_distance_map> Heat_method;
typedef CGAL::Heat_method_3::Heat_method_Eigen_traits_3::SparseMatrix SparseMatrix;



void source_set_tests(Heat_method hm, const Mesh& sm)
{
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  const std::set<vertex_descriptor>& source_copy = hm.get_sources();
  assert(*(source_copy.begin()) == source);;
  assert(*(hm.sources_begin()) == source);
  assert(hm.remove_source(source));
  assert((hm.get_sources()).empty());
  assert(hm.add_source(*(vertices(sm).first)));
  assert(*(hm.sources_begin()) == source);
  assert(hm.add_source(*(std::next(vertices(sm).first,3))));
  assert(source != *(std::next(vertices(sm).first,3)));
  assert(*(hm.sources_begin()) == source);
  assert(!(hm.add_source(source)));
  assert(hm.remove_source(*(std::next(vertices(sm).first,3))));
}

void cotan_matrix_test(const SparseMatrix& c)
{
  double sum = 0;
  for(int k = 0; k<c.outerSize(); ++k)
  {
    for(SparseMatrix::InnerIterator it(c,k); it; ++it)
    {
      sum +=it.value();
    }
  }
  std::cout<<"sum is: " << sum <<"\n";
  //Every row should sum up to 0, allow for slight error for large meshes
  assert(sum < 0.001);
}

void mass_matrix_test(const SparseMatrix& M)
{
  double sum = 0;
  for(int k = 0; k<M.outerSize(); ++k)
    {
      for(SparseMatrix::InnerIterator it(M,k); it; ++it)
      {
        sum +=it.value();
      }
    }
    //total Area matrix should be equal to the sum of all faces on the mesh
    //have to allow for the error because of rounding issues: Andreas might be able to help with this?
    //this will only work for the pyramid mesh
    assert((sum-1.866025)<=0.000005);
}

void check_for_zero(const Eigen::VectorXd& u)
{
  for(int c_i = 0; c_i<4; c_i++)
  {
      assert(u(c_i,0)<0.00001);
  }
}

void check_for_unit(const Eigen::MatrixXd& X, int dimension)
{
  for(int k = 0; k<dimension; k++)
  {
    double sum = CGAL::sqrt(X(k,0)*X(k,0) + X(k,1)*X(k,1) + X(k,2)*X(k,2));
    assert((sum-1)<0.00001);
  }
}

void check_no_update(const Mesh& sm, const Vertex_distance_map& original, const Vertex_distance_map& updated)
{
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm))
  {
    assert(get(original, vd) == get(updated,vd));
  }
}




int main()
{
  Mesh sm;
  Vertex_distance_map vertex_distance_map = get(Vertex_distance_tag(),sm);
  bool idf = false;

  std::ifstream in("../data/pyramid0.off");
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  //source set tests
  Heat_method hm(sm, vertex_distance_map, idf);
  source_set_tests(hm,sm);
  //cotan matrix tests
  const SparseMatrix& M = hm.mass_matrix();
  //std::cout<<"and M is: "<< Eigen::MatrixXd(M) << "\n";
  const SparseMatrix& c = hm.cotan_matrix();
  cotan_matrix_test(c);
  mass_matrix_test(M);

  double time_step = hm.time_step();
  double length_sum = hm.summation_of_edges();
  //there are 6 edges in pyramid
  double time_step_computed = (1./6)*length_sum;
  assert(time_step_computed*time_step_computed ==time_step);


  const SparseMatrix& K = hm.kronecker_delta();
  // AF: I commented the assert as I commented in build()
  assert(K.nonZeros()==1);
  Eigen::VectorXd solved_u = hm.solve_cotan_laplace(M,c,K,time_step,4);
  Eigen::VectorXd check_u = ((M+time_step*c)*solved_u)-K;
  //check_for_zero(check_u);
  Eigen::MatrixXd X = hm.compute_unit_gradient(solved_u);
  check_for_unit(X,3);

  const SparseMatrix& XD = hm.compute_divergence(X,4);

  Eigen::VectorXd solved_dist = hm.solve_phi(c, XD,4);

  Heat_method hm_idt(sm, vertex_distance_map, true);
  source_set_tests(hm_idt, sm);
  const SparseMatrix& M_idt = hm_idt.mass_matrix();
  const SparseMatrix& c_idt = hm_idt.cotan_matrix();
  cotan_matrix_test(c_idt);
  mass_matrix_test(M_idt);
  time_step = hm_idt.time_step();
  length_sum = hm_idt.summation_of_edges();
  time_step_computed = (1./6)*length_sum;

  const SparseMatrix& K_idt = hm_idt.kronecker_delta();
  solved_u = hm_idt.solve_cotan_laplace(M_idt,c_idt,K_idt,time_step, 4);
  check_u =((M_idt-time_step*c_idt) *solved_u)-K_idt;
  check_for_zero(check_u);
  X=hm_idt.compute_unit_gradient(solved_u);

  const SparseMatrix& XD_idt = hm_idt.compute_divergence(X,4);
  Eigen::VectorXd solved_idt_dist = hm_idt.solve_phi(c_idt, XD_idt,4);
  std::cout<< (solved_dist-solved_idt_dist) << "\n";

  Mesh sm2;
  Vertex_distance_map vertex_distance_map_2 = get(Vertex_distance_tag(),sm2);
  std::cout<<"bunny time\n";
  std::ifstream in2("../data/sphere.off");
  in2 >> sm2;
  if(!in2 || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  //  vertex_descriptor v9631 = CGAL::SM_Vertex_index(9631);
  //  std::cout<<"and vd is: "<< v9631 << "\n";
  Heat_method hm2(sm2, vertex_distance_map_2, true);
  source_set_tests(hm2, sm2);
  const SparseMatrix& M2 = hm2.mass_matrix();
  const SparseMatrix& c2 = hm2.cotan_matrix();
  cotan_matrix_test(c2);
//  hm2.add_source(v9631);
  std::cout<<"start of file disk distances\n";
  hm2.update();
  const Eigen::VectorXd& solved_dist_disk = hm2.distances();
  Eigen::VectorXd lib_geo_disk(43562,1);
  std::string line;

  std::ifstream in3("../data/sphere.dists");
  if(!in3) //Always test the file open.
   {
     std::cerr << "Problem loading the input data" << std::endl;
     return 1;
   }
   int i = 0;
   while (std::getline(in3, line))
   {
       lib_geo_disk(i,0) = std::stod(line);
       i++;
   }
   std::cout<<"AND ErRoR IS: "<< (lib_geo_disk-solved_dist_disk) << "\n";







  return 0;
}
