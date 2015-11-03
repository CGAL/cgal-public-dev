#include <CGAL/Lifting_kernel_d/sort_swap.h>
//#include <CGAL/Lifting_kernel_d/determinant.h>

// The boost implementation of hash tables appeared in version 1.36. If the
// installed version is older, we abort compilation.
#include <boost/version.hpp>
#if BOOST_VERSION < 103600
#error This needs Boost 1.36 or newer
#endif

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

#include <Eigen/Eigen>

#ifndef CGAL_HASH_TABLE_SIZE_LIMIT
// This value comes from empirical observations. When working with rational
// points, coming from doubles, this limits the process size to around 2Gb.
#define CGAL_HASH_TABLE_SIZE_LIMIT 7999999
#endif

namespace CGAL{

template <class _K,
          typename _P=typename _K::Point_d,
          typename _FT=typename _K::FT>
class Lifting_kernel_d:
public _K{
        private:
        typedef _K                                              Kernel_d;
        typedef _FT                                             FT;
        public:
        typedef _P                                              Point_d;
        private:
        typedef std::vector<const Point_d&>                     Index;
        typedef boost::unordered_map<Index,FT>                  Table;
        typedef typename Table::iterator                        It;

        public:
        struct Orientation_d;
};

template <class _K,typename _P,typename _FT>
struct Lifting_kernel_d<_K,_P,_FT>::Orientation_d{

    typedef _K                                              Kernel_d;
    typedef CGAL::Real_embeddable_traits<FT>                RET;
    typedef typename RET::Sgn                               Sgn;


    // This function calls the Orientation_d predicate of the base
    // kernel (i.e., it is a non-lifting predicate). Its interface is
    // described in the Kernel_d manual, it needs d+1 points.
    template <class PointInputIterator>
    CGAL::Orientation
    operator()(PointInputIterator first,PointInputIterator last)const{
            return Kernel_d().orientation_d_object()(first,last);
    }

    // This is the lifting operator
    // It needs d+2 points and a lifting coordinate for each
    template <class PointInputIterator,class LiftingInputIterator>
    CGAL::Orientation
    operator()(PointInputIterator pfirst,PointInputIterator plast,
               LiftingInputIterator lfirst,LiftingInputIterator llast)const{
        
        int num_of_points = static_cast<int>(std::distance(pfirst,plast));
        int dim = pfirst->dimension();
        
        // range contains d+2 points of dimension d
        CGAL_assertion_msg(dim + 2 == num_of_points,
          "Lifted Orientation_d: needs first->dimension() + 2 many points.");
        int num_of_lifts = static_cast<int>(std::distance(lfirst,llast));
        // range contains d+2 liftings one for each point 
        CGAL_assertion_msg(num_of_points == num_of_lifts,
          "Lifted Orientation_d: needs num_of_points as many as  num_of_lifts.");
        
        FT det = 0;
        for(int idx=0; idx<num_of_points; ++idx){
          FT minor;
          if(0){//Check if the matrix is precomputed in the hash table
            //TODO: implement hashtable
          }else{
            //Create eigen matrix and compute the determinant
            //TODO: instead of constructing d+2 can we only construct one?
            typedef Eigen::Matrix<FT,Eigen::Dynamic,Eigen::Dynamic> MT;
            int msize = num_of_points - 1;
            MT m(msize,msize);
            int j=0;
            for(PointInputIterator pit=pfirst; pit<plast; ++pit){
              if(static_cast<int>(std::distance(pfirst,pit))!=idx){
                int i=0;
                for(typename Point_d::Cartesian_const_iterator 
                    cit=pit->cartesian_begin(); cit<pit->cartesian_end(); ++cit){
                  m(i++,j) = *cit;
                }
                m(i,j++) = 1;
              }  
            }  
            //std::cout << m << std::endl << std::endl;
            minor = m.determinant();
            //std::cout << minor << std::endl;
          }
          minor *= *(lfirst+idx);
          //std::cout << minor << std::endl;
          //TODO: store minor in hash table and associate it with the set of indices 
          det += ((idx+1)+(num_of_points-1))%2==0 ? minor : -1 * minor;
        }
        // To be consistent with original kernel if the dimension of the lifted 
        // points ie. dim+1 is odd change sign
        det *= (dim+1)%2==0 ? 1 : -1;
        //std::cout << "determinant=" << det << std::endl;
        return Sgn()(det);
    }
};

} // namespace CGAL
