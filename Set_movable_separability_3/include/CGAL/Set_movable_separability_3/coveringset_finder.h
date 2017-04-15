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
#include <CGAL/intersections.h>

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
       * h1: a b c   => a*x+b*y+c*z>0
       * h2: A B C   => A*x+B*y+C*z>0
       * (pre: c>0)
       */
      template <typename Kernel>
      enum HemisphereInteractionState hemisphereInteraction(const typename Kernel::Plane_3 &h1, const typename Kernel::Plane_3 &h2)
      {
	if(h2.c()==0)
	  return HemisphereInteractionState_Intersect;
	if(h2.a()*h1.c() == h2.c()*h1.a()  &&  h2.b()*h1.c() == h2.c()*h1.b())
	  {
	    if(h2.c()>0)
	      return Hemisphere_equal;
	    return Hemisphere_antipodel;
	  }
	return HemisphereInteractionState_Intersect;
      }
      /*
       * This function projects h to a plane parallel to projTo in its negative part.
       * (pre: projTo has a positive z coefficient, both planes pass through the origin)
       *
       * main Idea: intersect hemisphere h and the plane projTo. project this plane to the xy-plane
       */
      template <typename Kernel>
      inline typename Kernel::Line_2 projectHemisphereToPlaneAndReturnItsComplementary(const typename Kernel::Plane_3 &h,const typename Kernel::Plane_3 &projTo)
      {
	/*
	 *
	 * plane parallel to projTo on its negative side:
	 *  	Ax+By+Cz=-1 (C>0) => Cz=-(Ax+By+1)
	 *
	 *	h: 	ax+by+cz>0	=> 	Cax+Cby+cCz>0
	 *	__________________________________________________
	 * (Ca - cA)x + (Cb - cB)y  -  c  >0
	 *
	 * but we want its Complementary ,
	 * so (cA - Ca)x + (cB - Cb)y  +  c
	 */
	typedef typename Kernel::Line_2 Line_2;
	return Line_2(h.c()*projTo.a()-h.a()*projTo.c(),
		      h.c()*projTo.b()-h.b()*projTo.c(),
		      h.c());

      }

      /*
       * input: hemispheres -  hemisphere of hemisphere in format of planes that its positive part is the hemisphere
       * (reminder: if the plane is given in this format: ax+by+cz+d=0 - its positive part is ax+by+cz+d>0.
       * 	 that planes ax+by+cz+d=0 and -ax-by-cz-d=0 are antipodel to each other)
       *
       * and out - empty array of size 6 at least
       *
       * output: a covering set in out (notice.. "a" covering set... not "the" covering set)
       * return value: size of this covering set, there might be one out direction that is not covered (and its antipodal direction)
       *
       * notice that it is assumed that the hemispheres are created by a polyhedron facets, if not, the code might fail.
       * for an instance, all planes must pass through the origin)
       */
      template <typename Kernel>
//	uint8_t findCoveringSet(typename Kernel::FT hemispheres[][3], unsigned int n, unsigned int * out,typename Kernel::Point_3 * outDirection)
      uint8_t findCoveringSet(const std::vector<typename Kernel::Plane_3> hemispheres, unsigned int * out,typename Kernel::Point_3 * outDirection)
			      //typename Kernel::Direction_3 * outDirection)
      {
	typedef typename Kernel::Line_2 Line_2;
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
	 * 5. find 3 hemispheres that covers -h
	 * 	& alg in the paper - basically, we project each hemisphere in hemispheres a half plane on the plane tangent to -h
	 * 	  and take the complementary half plane. Now we find 3 half-planes that their intersection is empty &
	 * 6. if 2 hemisphere are enough to cover -h (and one is not), this two must intersect on two points on c(h) and cover the rest of c(h)
	 *    just find two hemispheres that cover this points (we will actually return this two directions since a direction is easier to handle them
	 *    a possible top facet)
	 * 7. if you must use 3 hemispheres, they must also cover c(h).
	 */
	unsigned int hi=0;
	//jsut to make the code simpler, we will choose h such that its c is positive - one must exist
	for(;hi<hemispheres.size();++hi)
	  {
	    if(hemispheres[hi].c()>0)
	      break;
	  }
	bool removeHFromCoveringSet=false;
	std::vector<std::pair<Line_2,unsigned int> > halfplanes;
	halfplanes.reserve(hemispheres.size() -1);

	for(int i;i<hemispheres.size();++i)
	  {
	    if(i==hi)continue;

	    switch(hemisphereInteraction(hemispheres[hi],hemispheres[i]))
	    {
	      case HemisphereInteractionState_Intersect:
		halfplanes.push_back(
		    std::make_pair(
			projectHemisphereToPlaneAndReturnItsComplementary(hemispheres[i],hemispheres[hi]),
				   i));
		break;
	      case Hemisphere_equal:
		removeHFromCoveringSet=true;
		break;
	      case Hemisphere_antipodel:
		//ADD CODE
		break;
	    }
	  }
	uint8_t coveringSetSize = find3Witnesses(halfplanes,out);
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


	    typedef typename Kernel::Intersect_3 Intersect_3;
	    typedef typename Kernel::Line_3 Line_3;

	    Line_3 l = boost::get<Line_3>(* intersection(hemispheres[out[0]], hemispheres[out[1]]));
	    *outDirection= l.point(1);//get some arbitrary point in the intersection
	    if(outDirection->x()==0 && outDirection->y()==0 && outDirection->z()==0) //if it is the origin, take another point
	      *outDirection=l.point(2);
	    return coveringSetSize;
	  }


      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_COVERINGSET_FINDER_H_ */
