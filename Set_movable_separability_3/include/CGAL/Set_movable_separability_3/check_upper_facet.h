/*
 * check_upper_facet.h
 *
 *  Created on: Jul 8, 2017
 *      Author: shahar
 */

#ifndef INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_CHECK_UPPER_FACET_H_
#define INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_CHECK_UPPER_FACET_H_
/*
 * Unlike most of single_mold_translational_casting_3, this code is not based on
 * the paper "On the Separation of a Polyhedron from Its Single-Part Mold".
 * This code is based on chapter 4 in "computational geometry algorithms and applications"
 *
 * Given a polyhedron and a specific top facet, we are going to find whether or
 * not this top facet is valid. and if it is valid we are going to:
 *
 * 	 find a single pull-out direction for it in O(n) time.
 * OR
 * 	 find a all pull-out directions for it in O(n*logn) time.
 */
#include "Utils.h"
#include "PlaneProjector.h"
#include "lp_wrapper.h"

namespace CGAL {
  namespace Set_movable_separability_3 {
    namespace internal {
      template <typename Kernel,typename UserData>
         std::pair<bool,unsigned int> checkDirection(
             const std::vector< std::pair<typename Kernel::Direction_3, UserData>> &outerNormals,
						  typename Kernel::Direction_3 d)
         {
	int i;
	  for(i=0;i<outerNormals.size();++i)
	    {
	      if(outerNormals[i].first.dx()*d.dx()+
		  outerNormals[i].first.dy()*d.dy()+
		  outerNormals[i].first.dz()*d.dz()> 0)
		{


		  for(int j=i+1;j<outerNormals.size();++j)
		  	    {
		      if(outerNormals[j].first.dx()*d.dx()+
		    		  outerNormals[j].first.dy()*d.dy()+
		    		  outerNormals[j].first.dz()*d.dz()> 0)
		    		{

				  return  std::make_pair(false,UINT32_MAX);
		    		}
		  	    }
		  return  std::make_pair(true,i);
		}


	    }
	  std::cout<<"ERROR"<<std::endl;
	  return  std::make_pair(false,UINT32_MAX);//shouldn't get here with a valid poly
         }


      template <typename Kernel, typename UserData>
      std::pair<bool,typename Kernel::Direction_3> checkUpperFacet(
	  const   std::vector< std::pair<typename Kernel::Direction_3, UserData>> &outerNormals,
					unsigned int iTopFacet)
      {
	typedef typename Kernel::Direction_3 Direction_3;

	typedef typename Kernel::Line_2 Line_2;
	PlaneProjector<Kernel,false> proj(outerNormals[iTopFacet].first);

	  int i;
	  std::vector<Line_2> halfplanes;
	  halfplanes.reserve(outerNormals.size()-1);
	  for(i = 0;i<outerNormals.size();++i)
	    {
	      if(unlikely(i==iTopFacet))
		continue;

	      std::pair<enum LineState,Line_2> line = proj.projectHemisphereToPlaneAndReturnItsComplementary(outerNormals[i].first,false);
	      if(likely(line.first==LINESTATE_LINE))
		{
		  halfplanes.push_back(line.second);
		}
	      if(line.first==LINESTATE_SAME_NORAML)
		return  std::make_pair(false,Direction_3());


	      if(line.first==LINESTATE_ANTIPODEL_NORMAL)
		continue;

	    }
	       std::pair<bool, typename Kernel::Point_2> pulloutDirection =  findPoint<Kernel>(halfplanes);
	       if(pulloutDirection.first)
		 {
			return  std::make_pair(true,proj.point_2ToDirection_3(pulloutDirection.second));

		 }
	       else
		 {
			return  std::make_pair(false,proj.point_2ToDirection_3(pulloutDirection.second));

		 }

      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_CHECK_UPPER_FACET_H_ */
