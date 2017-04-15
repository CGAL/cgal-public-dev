/*
 * coveringsetfinder.h
 *
 *  Created on: Apr 14, 2017
 *      Author: shahar
 */

#ifndef INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_COVERINGSET_FINDER_H_
#define INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_COVERINGSET_FINDER_H_
/*
 * This function goal is to find a "covering-set" of size 6 or less for a polyhedron.
 * (full definition and motivation can be found in our paper:
 * "On the Separation of a Polyhedron from Its Single-Part Mold" (version of 2017)
 *
 * Semi-definition:
 * * direction: a point of the unit-sphere (S^2) that represent the direction from the origin to it.
 * * The directions of a facet F in a polyhedron: all the directions that if you move epsilon along them from a point in the
 *   middle of the facet, you get leave the polyhedron. Notice that this is an open hemisphere (in the paper it called h_bar).
 * * covering-set: a set of hemisphere that the union of their directions, is the entire unit sphere.
 *   (a proof of the existence of such set of size 6 - in the paper)
 */
#include "lp_wrapper.h"
#include <vector>
namespace CGAL {
  namespace Set_movable_separability_3 {
    namespace internal {

      enum HemisphereInteractionState
      {
	HemisphereInteractionState_Intersect,
	Hemisphere_equal,
	Hemisphere_antipodel
      };
      /*
       * a1 b1 c1   => a1*x+b1*y+c1*z<0
       * a2 b2 c2   => a2*x+b2*y+c2*z<0
       */
      template <typename FT>
      enum HemisphereInteractionState hemisphereInteraction(FT *h1,FT* h2)
      {
	if(h1[0]==0)
	  {
	    if(h2[0]!=0)
	      return HemisphereInteractionState_Intersect;
	    if(h1[1]*h2[2]!=h1[2]*h2[1])
	      return HemisphereInteractionState_Intersect;
	    if(h1[1]==0)
	      {
		if( (h2[2]>0) == (h2[2]>0) )
		  return Hemisphere_equal;
		return Hemisphere_antipodel;
	      }
	    if( (h2[1]>0) == (h2[1]>0) )
	      return Hemisphere_equal;
	    return Hemisphere_antipodel;

	  }
	else
	  {
	    if(h2[0]==0)
	      return HemisphereInteractionState_Intersect;
	    if(h1[1]*h2[0]!=h1[0]*h2[1] ||  h1[2]*h2[0]!=h1[0]*h2[2])
	      return HemisphereInteractionState_Intersect;

	    if( (h2[0]>0) == (h2[0]>0) )
	      return Hemisphere_equal;
	    return Hemisphere_antipodel;
	  }
      }

      /*
       * input: hemispheres -  array of hemisphere in format:
       * a1 b1 c1   => a1*x+b1*y+c1*z<0
       * a2 b2 c2   => a2*x+b2*y+c2*z<0
       * a3 b3 c3   => a3*x+b3*y+c3*z<0
       * ....
       * an bn cn   => an*x+bn*y+cn*z<0
       *
       * and out - empty array of size 6 at least
       *
       * output: a covering set in out (notice.. "a" covering set... not "the" covering set)
       * return value: size of this covering set, there might be one out direction that is not covered (and its antipodal direction)
       *
       * notice that it is assumed that the hemispheres are created by a polyhedron facets, if not, the code might fail.
       */
      template <typename Kernel>
	uint8_t findCoveringSet(typename Kernel::FT hemispheres[][3], unsigned int n, unsigned int * out,typename Kernel::Point_3 * outDirection)
			      //typename Kernel::Direction_3 * outDirection)
      {
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Point_3 Point_3;
	/*
	 * algorithm: this alg is not written in the current version of our paper, so I add it in here.
	 * (pre: the union of all hemispheres is S^2 - will occur in any polyhedron)
	 *
	 * symbols:
	 * h - open hemisphere
	 * -h - the parallel open hemisphere
	 * c(h) - the great circle of h
	 * (notice that the union of h, -h and c(h) is S^2 for any h)
	 *
	 * 1. covering-set <- {}
	 * 2. h <- choose some arbitrary hemisphere of hemispheres
	 * 3. covering-set.add(h)
	 *	& -h and c(h) &
	 * 4. if some hemisphere in hemispheres is -h, we just need to cover c(h) with 3 hemispheres and finish (easy, I won't explain it)
	 * 5. find 3 hemispheres that covers -h (alg in the paper)
	 * 6. if 2 hemisphere are enough to cover -h (and one is not), this two must intersect on two points on c(h) and cover the rest of c(h)
	 *    just find two hemispheres that cover this points (we will actually return this two directions since a direction is easier to handle them
	 *    a possible top facet)
	 * 7. if you must use 3 hemispheres, they must also cover c(h).
	 */
	unsigned int hi=0;
	//jsut to make the code easier, we will chose h such that its c is positive - one must exist
	for(;hi<n;++hi)
	  {
	    if(hemispheres[hi][2]>0)
	      break;
	  }
	int halfPlaneCount=0;
	bool removeHFromCoveringSet=false;
	FT halfplanes[n-1][3];

	for(int i;i<n;++i)
	  {
	    if(i==hi)continue;

	    switch(hemisphereInteraction(hemispheres[i][2],hemispheres[hi][2]))
	    {
	      case HemisphereInteractionState_Intersect:
		halfplanes[halfPlaneCount][0]= hemispheres[hi][2]* hemispheres[i][0] -hemispheres[hi][0]* hemispheres[i][2];
		halfplanes[halfPlaneCount][1]= hemispheres[hi][2]* hemispheres[i][1] -hemispheres[hi][1]* hemispheres[i][2];
		halfplanes[halfPlaneCount][2]= -hemispheres[i][2];
		halfPlaneCount++;
		break;
	      case Hemisphere_equal:
		removeHFromCoveringSet=true;
		break;
	      case Hemisphere_antipodel:
		//ADD CODE
		break;
	    }
	  }
	uint8_t coveringSetSize = find3Witnesses(halfplanes,halfPlaneCount,out);
	if(coveringSetSize==3)
	  {
	    if(!removeHFromCoveringSet)
	      out[coveringSetSize++]=hi;
	    return coveringSetSize;
	  }
	else //coveringSetSize==2
	  {
	    if(!removeHFromCoveringSet)
	      out[coveringSetSize++]=hi;
	  //  CGAL::intersection()
	    return coveringSetSize;
	  }


      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_COVERINGSET_FINDER_H_ */
