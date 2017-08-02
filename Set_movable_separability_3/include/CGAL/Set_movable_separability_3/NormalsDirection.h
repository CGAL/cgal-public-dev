/*
 * NormalsDirection.h
 *
 *  Created on: Jun 14, 2017
 *      Author: shahar
 */

#ifndef INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_NORMALSDIRECTION_H_
#define INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_NORMALSDIRECTION_H_
#include <vector>

/*
 * This function's goal is to figure out for each of the polyhedron's faces to the faces, the normal that
 * is oriented toward the outer direction
 */


#include "Utils.h"

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//
//typedef CGAL::Polyhedron_3<Kernel>                          Poly;

namespace CGAL {
  namespace Set_movable_separability_3 {
    namespace internal {
#include<vector>
#include <iostream>


      template<typename Poly,typename Kernel>
      std::vector< typename Kernel::Direction_3> findDirections(const Poly& poly)
      {
	typedef typename Poly::Halfedge_const_handle  Halfedge_const_handle;
        typedef typename Kernel::Direction_3 Direction_3;
        typedef typename Kernel::FT FT;
        std::vector<Direction_3> nd;
        nd.reserve(poly.size_of_facets());
        for(typename Poly::Facet_const_iterator fit= poly.facets_begin();fit != poly.facets_end();++fit)
        {

            Halfedge_const_handle heCurrent=fit->halfedge();
            Halfedge_const_handle heBeforeMax=heCurrent;

            heCurrent= heCurrent->next();
            FT max=heCurrent->vertex()->point().z();
            bool allEqual=true;

            Halfedge_const_handle heStart=heCurrent;
            Halfedge_const_handle heBeforeCurrent=heCurrent;
            heCurrent= heCurrent->next();

            int countDbg=1;
            int countEqual=1;

            do
            {
        	if(heCurrent->vertex()->point().z()>max)
        	{
        	    max=heCurrent->vertex()->point().z();
        	    heBeforeMax=heBeforeCurrent;
        	    allEqual=false;
        	}
        	else if(unlikely(allEqual))
        	{
        	    countEqual++;
        	    allEqual = heCurrent->vertex()->point().z()==max;
        	    if(allEqual && countEqual==3)
        	      break;
        	}
        	heCurrent= heCurrent->next();
            }while (heCurrent!=heStart);
            if(unlikely(allEqual))
	    {
		 heCurrent=fit->halfedge();
		 heCurrent= heCurrent->next();
		 max=heCurrent->vertex()->point().x();
		 heBeforeCurrent=heCurrent;
		 heCurrent= heCurrent->next();

		  do
		  {
		      if(heCurrent->vertex()->point().x()>max)
		      {
			  max=heCurrent->vertex()->point().x();
			  heBeforeMax=heBeforeCurrent;
		      }
		      heCurrent= heCurrent->next();
		  }while (heCurrent!=heStart);
	    }
typedef typename Poly::Point Point;
    	Point a= heBeforeMax->vertex()->point();
    	Point  b=heBeforeMax->next()->vertex()->point();
    	Point  c= heBeforeMax->next()->next()->vertex()->point();

    	Direction_3 nor(CGAL::normal(a,b,c));
//    	std::cout<<"["<<a<< "|"<<b<< "|"<<c<< "] -> "<<nor<<std::endl;
    	nd.push_back(nor);
        }
        return nd;
      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL


#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_NORMALSDIRECTION_H_ */
