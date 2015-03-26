#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <omp.h>



int main(int argc,char** argv)
{
  omp_set_num_threads( omp_get_max_threads() );
  
  CGAL::Gmpq n1(1);
  CGAL::Gmpq cumul(0);
  
  int loop_nb=1000;
  if (argc!=1)
    loop_nb=atoi(argv[1]);
  
  #pragma omp parallel for
  for(int x=0; x < loop_nb; x++)
  { 
    CGAL::Gmpq n2(n1);
    
    #pragma omp critical
    cumul+=n2;
  }

  std::cout << cumul << std::endl;
}


