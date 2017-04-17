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
 *
 * Motivation: This helps as find all the top facets. We know that the hemisphere of a top facet must be the only hemisphere
 * 	       that covers some cell on S^2. By that we know that if we will find a covering set - all the top-facets hemispheres
 * 	       must be in this set.

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
       * This function projects h to a line parallel to projTo in its negative part
       * (pre: projTo has a positive y coefficient, both line pass through the origin and do not merge)
       *
       * main Idea: intersect halfcircle h and the line projTo. project this plane to the x-axe
       *
       * return pair (x threshold value, true for bigger than this x - false for smaller)
       */
      template <typename Kernel>
      inline std::pair<typename Kernel::FT,bool> projectHalfcircleToline(const typename Kernel::Line_2 &h,const typename Kernel::Line_2 &projTo)
      {
	/*
	 *
	 * line parallel to projTo on its negative side:
	 *  	Ax+By=-1 (B>0) => By=-(Ax+1)
	 *
	 *	h: 	ax+by>0	=> 	Bax+Bby>0
	 *	__________________________________________________
	 * (Ba - bA)x + -b  >0
	 * (Ba - bA)x > b
	 * threshold is b/(Ba - bA)
	 * if (Ba - bA)>0 we want bigger values
	 * if (Ba - bA)<0 we want smaller values
	 * notice that (Ba - bA)!=0 since the lines are different
	 */
	typedef typename Kernel::FT FT;
	FT Ba_bA = projTo.b()*h.a()- projTo.a()*h.b();
	return std::make_pair(h.b()/Ba_bA,Ba_bA>0);
      }

      /*
       * This function projects h to a projTo OR plane parallel to projTo in its negative part (by parameter)
       * (pre: projTo has a positive z coefficient, both planes pass through the origin and do not merge)
       *
       * main Idea: intersect hemisphere h and the plane projTo. project this plane to the xy-plane
       */
      template <typename Kernel>
      inline typename Kernel::Line_2 projectHemisphereToPlaneAndReturnItsComplementary(const typename Kernel::Plane_3 &h,const typename Kernel::Plane_3 &projTo,bool projToNegative)
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
		      ( projTo?(h.c()):0)
		      );

      }
      /*
       * pretty much the same as the function below (findCoveringSet)
       * just in 2D.

       * input:
       * * halfcircle -  vector of halfcircles in format of Line_2 that it's positive part is the halfcircle
       *   (reminder: if the line is given in this format: ax+by+c=0 - its positive part is ax+by+c>0.
       * 	 Those lines: ax+by+cz+d=0 and -ax-by-cz-d=0 are antipodel to each other
       * 	 In our case c is always 0 )
       * 	 Each halfcircle is paired with its name as int.
       *
       * * out - empty array of size 3 unsigned int at least
       *
       * (outDirection - is not given in here since we don't have the original hemispheres thus we can't compute
       * 		 the single direction we didn't cover. if our return value is 2, the caller function need to compute this direction)
       *
       * output: a covering set names (the paired int) in out (notice.. "a" covering set... not "the" covering set)
       * return value: size of this covering set,there might be one out direction that is not covered (and its antipodal direction).
       * 		in this case isOutDirection is set to true.
       *
       * This is not actually true... since the covering-set in only a means to an end (read the motivation above), we will not return
       * halfcircle that we noticed that appears twice - so the return value might not be a covering-set (that doesn't mean that we checked it for any halfplane returned).
       *
       * notice that it is assumed that the halfcircle are created by hemispherescreated by a polyhedron facets, if not, the code might fail.
       * for an instance, all lines must pass through the origin.

       */
      template <typename Kernel>
          uint8_t findCoveringSet2D(const std::vector<std::pair<typename Kernel::Line_2,unsigned int>>  halfcircle,
				    unsigned int * out,bool* isOutDirection)
      {
	typedef typename Kernel::Line_2 Line_2;
	typedef typename Kernel::FT FT;

	unsigned int hi=0;

	//just to make the code simpler, we will choose h such that its b is positive - one must exist
	//otherwise (0,1) wasn't in any hemisphere
	for(;hi<halfcircle.size();++hi)
	  {
	    if(halfcircle[hi].first.b()>0)
	      break;
	  }
	  bool removeHFromCoveringSet=false;

	  bool minRayToInfHasValue=false;
	  FT minRayToInf;
	  bool maxRayFromMinusInfHasValue=false;
	  FT maxRayFromMinusInf;
	  for(int i;i<halfcircle.size();++i)
	    {
	      if(i==hi)continue;

	      switch(halfPlaneInteraction(halfcircle[hi].first,halfcircle[i].first))
	      {
		case HalfplaneInteractionState_LinesIntersect:
		  std::pair<FT,bool> ray = projectHalfcircleToline(halfcircle[i],halfcircle[hi]);
		   if(ray.second)
		     //bigger than
		     {
		       if(!minRayToInfHasValue)
			 {
			   minRayToInfHasValue=true;
			   minRayToInf=ray.first;
			   out[0]=halfcircle[i].second;
			 }
		       else
			 {
			     if(minRayToInf>ray.first)
			       {
				 minRayToInf=ray.first;
				 out[0]=halfcircle[i].second;
			       }
			 }
		     }
		   else
		     //smaller  than
		     {
		       if(!maxRayFromMinusInfHasValue)
			 {
			   maxRayFromMinusInfHasValue=true;
			   maxRayFromMinusInf=ray.first;
			   out[1]=halfcircle[i].second;
			 }
		       else
			 {
			   if(maxRayFromMinusInf<ray.first)
			     {
			       maxRayFromMinusInf=ray.first;
			       out[1]=halfcircle[i].second;
			     }
			 }
		     }
		   break;
		case HalfplaneInteractionState_LinesParallel_FisrtContainsSecond:case HalfplaneInteractionState_LinesParallel_SecondContainsFirst:
		  //we know that since all line goes throw the origin, both means equal
		  removeHFromCoveringSet=true;
		  break;
		case HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty://antipodel
		  out[0]=halfcircle[i].second;
		  *isOutDirection=true;
		  if(!removeHFromCoveringSet)
		    {
		      out[1]=halfcircle[hi].second;
		      return 2;
		    }
		  return 1;
		  break;
	      }
	    }
	  //assert minRayToInf<maxRayFromMinusInf
	  *isOutDirection=false;
	  if(!removeHFromCoveringSet)
	  {
	    out[2]=halfcircle[hi].second;
	    return 3;
	  }
	return 2;


      }
      /*
       * input:
       * * hemispheres -  vector of hemisphere in format of planes that it's positive part is the hemisphere
       *   (reminder: if the plane is given in this format: ax+by+cz+d=0 - its positive part is ax+by+cz+d>0.
       * 	 Those planes" ax+by+cz+d=0 and -ax-by-cz-d=0 are antipodel to each other
       * 	 In our case d is always 0)
       *
       * * out - empty array of size 6 unsigned int at least
       *
       * *outDirection - pointer to an empty Point_3
       *
       * output: a covering set in out (notice.. "a" covering set... not "the" covering set)
       * return value: size of this covering set,there might be one out direction that is not covered (and its antipodal direction)/.
       * if such direction exists, it is returned in outDirection. if not - outDirection is the origin.
       * This is not actually true... since the covering-set in only a means to an end (read the motivation above), we will not return
       * hemispheres that we noticed that appears twice - so the return value might not be a covering-set (that doesn't mean that we checked it for any hemisphere returned).
       *
       * notice that it is assumed that the hemispheres are created by a polyhedron facets, if not, the code might fail.
       * for an instance, all planes must pass through the origin.
       */
      template <typename Kernel>
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
	 * 4. if some hemisphere in hemispheres is -h, we just need to cover c(h) with 3 hemispheres and finish (Same idea - one less dimension)
	 * 5. find 3 hemispheres that covers -h
	 * 	& alg in the paper - basically, we project each hemisphere in hemispheres a half plane on the plane tangent to -h
	 * 	  and take the complementary half plane. Now we find 3 half-planes that their intersection is empty &
	 * 6. if 2 hemisphere are enough to cover -h (and one is not), this two must intersect on two points on c(h) and cover the rest of c(h)
	 *    just find two hemispheres that cover this points (we will actually return this two directions since a direction is easier to handle them
	 *    a possible top facet)
	 * 7. if you must use 3 hemispheres, they must also cover c(h).
	 */
	unsigned int hi=0;
	//just to make the code simpler, we will choose h such that its c is positive - one must exist
	//otherwise (0,0,1) wasn't in any hemisphere
	for(;hi<hemispheres.size();++hi)
	  {
	    if(hemispheres[hi].c()>0)
	      break;
	  }
	bool removeHFromCoveringSet=false;
	bool hBarFound=false;
	bool removeHBarFromCoveringSet;
	bool hBarIndex;

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
			      projectHemisphereToPlaneAndReturnItsComplementary(hemispheres[i],hemispheres[hi],
			      !hBarFound ),//if HbarFound we want to find a covering set for c(h), else for -h
			i));
		break;
	      case Hemisphere_equal:
		removeHFromCoveringSet=true;
		break;
	      case Hemisphere_antipodel:
		removeHBarFromCoveringSet = hBarFound;//if there is more than one hBar, this cant be a top facet
		hBarFound=true;
		hBarIndex=i;
		for(int j = 0; j!=halfplanes.size();j++)
		  {
		    halfplanes[j].first=Line_2(halfplanes[j].first.a(),halfplanes[j].first.b(),0);
		  }
		break;
	    }
	  }
	    typedef typename Kernel::Intersect_3 Intersect_3;
	    typedef typename Kernel::Line_3 Line_3;
	if(hBarFound)
	  {
	    bool isOutDirection=false;
	    uint8_t coveringSetSize =findCoveringSet2D(halfplanes,out,outDirection,&isOutDirection);
	    if(isOutDirection)
	      //this is only possible if the two lines chosen are the same line to opposite sides.
	      //in this case, the direction of this line isn't covered
	      {
		 Line_3 l = boost::get<Line_3>(* intersection(hemispheres[out[0]], hemispheres[hi]));
		*outDirection= l.point(1);//get some arbitrary point in the intersection
		if(outDirection->x()==0 && outDirection->y()==0 && outDirection->z()==0) //if it is the origin, take another point
		  *outDirection=l.point(2);
	      }
	    if(!removeHFromCoveringSet)
	      out[coveringSetSize++]=hi;
	    if(!removeHBarFromCoveringSet)
	      out[coveringSetSize++]=hBarIndex;
	    return coveringSetSize;

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
