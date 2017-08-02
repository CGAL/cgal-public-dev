/*
 * Utils.h
 *
 *  Created on: Jun 14, 2017
 *      Author: shahar
 */

#ifndef INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_UTILS_H_
#define INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_UTILS_H_


namespace CGAL {
  namespace Set_movable_separability_3 {
    namespace internal {
      #define likely(x)       __builtin_expect((x),1)
      #define unlikely(x)     __builtin_expect((x),0)
      typedef bool NoramlDirection; //for normal (a,b,c) true is a>0 || (a==0 && b>0) || (a==0 && b==0 && c>0)
/*
      template<typename Kernel>
      Kernel::Plane_3 threePointToDirectedPlane(typename Kernel::Point_3 firstInFacet,
						typename Kernel::Point_3 secondInFacet,
						typename Kernel::Point_3 thirdInFacet,
						NoramlDirection facetDirection)
      {

      }*/


    }
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL


#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_UTILS_H_ */
