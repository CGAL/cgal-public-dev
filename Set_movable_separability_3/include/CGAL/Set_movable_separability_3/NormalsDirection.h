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
#include<map>
#include <iostream>
      /*
       * output:
       * vector of outer normal , facet (facet is represented by a list of facet for the case that this facet is splited in to coplanar facets)
       */
      template<typename Poly,typename Kernel>
    	void innerMergeCoPlanarFacets(const Poly& poly,
    			       typename Poly::Facet_const_iterator fit,
    			       typename Kernel::Direction_3 normal,
    			        std::map<typename Poly::Facet_const_iterator,typename Kernel::Direction_3>&directions,
    				std::vector<typename Poly::Facet_const_iterator> & outPutFacetList
    			      )
          {
    	typedef typename Poly::Halfedge_const_handle  Halfedge_const_handle;
    	outPutFacetList.push_back(fit);
    	Halfedge_const_handle heCurrent=fit->halfedge();
    	Halfedge_const_handle heStart=heCurrent;


    	do
    	  {
    	    typename Poly::Facet_const_iterator neighbor = heCurrent->opposite()->facet();
    	   typename  std::map<typename Poly::Facet_const_iterator,typename Kernel::Direction_3>::iterator dit= directions.find(neighbor);

    	    if(dit != directions.end() && ((dit->second)==normal))
    	      {
    		directions.erase(dit);

    		innerMergeCoPlanarFacets<Poly,Kernel>(poly,neighbor,normal,directions,outPutFacetList);
    	      }
    	    heCurrent= heCurrent->next();
    	  }while (heCurrent!=heStart);


          }
      template<typename Poly,typename Kernel>
      std::vector< std::pair<typename Kernel::Direction_3, std::vector<typename Poly::Facet_const_iterator> > >
      mergeCoPlanarFacets(const Poly& poly, std::map<typename Poly::Facet_const_iterator,typename Kernel::Direction_3> &directions)
      {
//	      std::cout<<directions.size()<<" VS ";
	      std::vector< std::pair<typename Kernel::Direction_3, std::vector<typename Poly::Facet_const_iterator>>> normals;
	      for(typename Poly::Facet_const_iterator fit= poly.facets_begin();fit != poly.facets_end();++fit)
		{
		  typename std::map<typename Poly::Facet_const_iterator,typename Kernel::Direction_3>::iterator dit= directions.find(fit);

		  if(dit!=directions.end())
		    {
		      std::pair<typename Kernel::Direction_3, std::vector<typename Poly::Facet_const_iterator>> normal;
		      normal.first=dit->second;
		      directions.erase(dit);
		      innerMergeCoPlanarFacets<Poly,Kernel>(poly,fit,dit->second,directions, normal.second);
		      normals.push_back(normal);
		    }
		}
//	      std::cout<<normals.size()<<std::endl;
	      return normals;
      }

      template<typename Poly,typename Kernel>
      std::vector< std::pair<typename Kernel::Direction_3, std::vector<typename Poly::Facet_const_iterator> > > findDirections(const Poly& poly)
      {
	typedef typename Poly::Halfedge_const_handle  Halfedge_const_handle;
        typedef typename Kernel::Direction_3 Direction_3;
        typedef typename Kernel::FT FT;
        std::map<typename Poly::Facet_const_iterator,Direction_3> nd;
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
        	heBeforeCurrent= heCurrent;
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
		      heBeforeCurrent= heCurrent;
		      heCurrent= heCurrent->next();
		  }while (heCurrent!=heStart);
	    }
typedef typename Poly::Point Point;
	  Point a= heBeforeMax->vertex()->point();
	  Point  b=heBeforeMax->next()->vertex()->point();
    	  Point  c= heBeforeMax->next()->next()->vertex()->point();

    	  Direction_3 nor(CGAL::normal(a,b,c));
    	  nd[fit]=nor;
        }


        return mergeCoPlanarFacets<Poly,Kernel>(poly,nd);
      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL


#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_NORMALSDIRECTION_H_ */
