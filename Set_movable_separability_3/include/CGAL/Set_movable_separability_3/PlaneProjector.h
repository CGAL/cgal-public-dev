/*
 * PlaneProjector.h
 *
 *  Created on: Jul 25, 2017
 *      Author: shahar
 */

#ifndef INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_PLANEPROJECTOR_H_
#define INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_PLANEPROJECTOR_H_
#include "Utils.h"

namespace CGAL {
  namespace Set_movable_separability_3 {
    namespace internal {
	enum LineState{
	  LINESTATE_LINE,
	  LINESTATE_ALL_PLANE,
	  LINESTATE_EMPTY
	};
      template <typename Kernel>
      class PlaneProjectorFather
      {
      public:

	virtual ~PlaneProjectorFather(){}
	//returns false is this facet is not a top facet for sure
	virtual std::pair<enum LineState,typename Kernel::Line_2> projectHemisphereToPlaneAndReturnItsComplementary(typename Kernel::Direction_3 normal,bool projToPlaneThroughOrigin)=0;
	virtual typename Kernel::Direction_3 point_2ToDirection_3(typename Kernel::Point_2 p)=0;



      };
      template <typename Kernel, int TopFacetCSign,int TopFacetBSignIfCIsZero,int TopFacetASignIfBAndCAreZero,bool isOnNegetive>
      	class InnerPlaneProjector:public PlaneProjectorFather<Kernel>
      	{
      	  typename Kernel::Direction_3 m_normalOfMainPlane;//(A,B,C)
  	virtual ~InnerPlaneProjector(){}
      	public:

      	  InnerPlaneProjector(typename Kernel::Direction_3 normalOfMainPlane):m_normalOfMainPlane(normalOfMainPlane)
      	  {
      	  }
      	virtual typename Kernel::Direction_3 point_2ToDirection_3(typename Kernel::Point_2 p)
      	{
      	  /*
      	   * We want to find a point (/direction) (p.x,p.y,z) such that
      	   * A*p.x+B*p.y+C*z=+-1 (C!=0)
      	   * z=(-A*p.x-B*p.y+-1)/C
      	   */
      	  typedef typename Kernel::Direction_3 Direction_3;
	    if(TopFacetCSign!=0)//compile time if! how cool
    	      {
    		return  Direction_3(
    		    p.x(),
		    p.y(),
		    (-m_normalOfMainPlane.dx()*p.x()
		    -m_normalOfMainPlane.dy()*p.y()
		    +(isOnNegetive?-1:1))/
		    m_normalOfMainPlane.dz()
    		);
    	      }
    	    else if(TopFacetBSignIfCIsZero!=0)//&&(TopFacetCSign==0)//compile time if! how cool
    	      {
    		/*
    		 *  We want to find a point (/direction) (p.x,Y,p.y) such that
    		 * A*p.x+B*Y=+-1 (B!=0)
    		 * Y=(-A*p.x+-1)/B
    		 */

    		return  Direction_3(
    	    		    p.x(),

    			    (-m_normalOfMainPlane.dx()*p.x()
    			    +(isOnNegetive?-1:1))/
    			    m_normalOfMainPlane.dy(),
			    p.y()
    	    		);

    	      }
    	    else// if(TopFacetASignIfBAndCAreZero!=0) &&(TopFacetBSignIfCIsZero==0)&&(TopFacetCSign==0)//compile time if! how cool
    	      {
    		/*
    		 * We want to find a point (/direction) (X,p.x,p.y) such that
    		 * A*X=+-1 (B!=0)
    		 * X=+-1/A
    		 *
    		 * plane parallel to projTo on its negative side:
    		 *  	Ax=-1 (A>0)
    		 *
    		 *	h: 	ax+by+cz>0	=> 	aAx+Aby+Acz>0
    		 *	__________________________________________________
    		 * Aby+Acz -a >0
    		 *
    		 * but we want its Complementary,
    		 * -Aby-Acz +a >0
    		 */

    		return  Direction_3(
    			    (isOnNegetive?-1:1)/
    			    m_normalOfMainPlane.dx(),
    	    		    p.x(),
			    p.y()
    	    		);

    	      }


      	}

      	  //returns false is this facet is not a pot facet for sure
      	  virtual std::pair<enum LineState,typename Kernel::Line_2> projectHemisphereToPlaneAndReturnItsComplementary(typename Kernel::Direction_3 normal,bool projToPlaneThroughOrigin)
      	  {
      	    /* normal (a,b,c)
      	     *
      	     * plane parallel to projTo on its positive/negative side:
      	     *  	Ax+By+Cz=+-1 (C>0) => Cz=-(Ax+By+-1)
      	     *
      	     *	h: 	ax+by+cz>0	=> 	Cax+Cby+cCz>0
      	     *	__________________________________________________
      	     * (Ca - cA)x + (Cb - cB)y  +-  c  >0
      	     *
      	     * but we want its Complementary ,
      	     * so (cA - Ca)x + (cB - Cb)y  -+  c >0
      	     */
      	    typedef typename Kernel::Line_2 Line_2;
      	    if(unlikely(normal==m_normalOfMainPlane))
      	      return std::make_pair(LINESTATE_EMPTY,Line_2());
      	    if(unlikely(normal==-m_normalOfMainPlane))
      	      return std::make_pair(LINESTATE_ALL_PLANE,Line_2());
      	    if(TopFacetCSign>0)//compile time if! how cool
      	      {
      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(normal.dz()*m_normalOfMainPlane.dx()-normal.dx()*m_normalOfMainPlane.dz(),
      			normal.dz()*m_normalOfMainPlane.dy()-normal.dy()*m_normalOfMainPlane.dz(),
			projToPlaneThroughOrigin?0:(isOnNegetive?normal.dz():-normal.dz())
      		    ))
      		);
      	      }
      	    else if(TopFacetCSign<0)//compile time if! how cool
      	      { //same as the main comment just with minus (reverse sign)
      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(normal.dx()*m_normalOfMainPlane.dz()-normal.dz()*m_normalOfMainPlane.dx(),
      			normal.dy()*m_normalOfMainPlane.dz()-normal.dz()*m_normalOfMainPlane.dy(),
			projToPlaneThroughOrigin?0:(isOnNegetive?-normal.dz():normal.dz())
      		    ))
      		);
      	      }
      	    else if(TopFacetBSignIfCIsZero>0)//&&(TopFacetCSign==0)//compile time if! how cool
      	      {
      		/*
      		 * plane parallel to projTo on its positive/negative side:
      		 *  	Ax+By=+-1 (B>0) => By=-(Ax-+1)
      		 *
      		 *	h: 	ax+by+cz>0	=> 	Bax+bBy+Bcz>0
      		 *	__________________________________________________
      		 * (Ba - bA)x + Bcz +- b > 0
      		 *
      		 * but we want its Complementary,
      		 * (bA - Ba)x - Bcz -+ b > 0
      		 */


      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(normal.dy()*m_normalOfMainPlane.dx()-normal.dx()*m_normalOfMainPlane.dy(),
      			-normal.dz()*m_normalOfMainPlane.dy(),
			projToPlaneThroughOrigin?0:(isOnNegetive?normal.dy():-normal.dy())
      		    ))
      		);

      	      }
      	    else if(TopFacetBSignIfCIsZero<0)//&&(TopFacetCSign==0)//compile time if! how cool
      	      { //same as the last comment just with minus (reverse sign)
      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(normal.dx()*m_normalOfMainPlane.dy()-normal.dy()*m_normalOfMainPlane.dx(),
      			normal.dz()*m_normalOfMainPlane.dy(),
			projToPlaneThroughOrigin?0:(isOnNegetive?-normal.dy():normal.dy())))
      		);
      	      }
      	    else if(TopFacetASignIfBAndCAreZero>0) //&&(TopFacetBSignIfCIsZero==0)&&(TopFacetCSign==0)//compile time if! how cool
      	      {
      		/*
      		 * plane parallel to projTo on its negative side:
      		 *  	Ax=-1 (A>0)
      		 *
      		 *	h: 	ax+by+cz>0	=> 	aAx+Aby+Acz>0
      		 *	__________________________________________________
      		 * Aby+Acz -a >0
      		 *
      		 * but we want its Complementary,
      		 * -Aby-Acz +a >0
      		 */
      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(-normal.dy()*m_normalOfMainPlane.dx(),
      			-normal.dz()*m_normalOfMainPlane.dx(),
			projToPlaneThroughOrigin?0:(isOnNegetive?normal.dx():-normal.dx())
			))
      		);

      	      }
      	    else // if(TopFacetASignIfBAndCAreZero<0 &&(TopFacetBSignIfCIsZero==0)&&(TopFacetCSign==0))//compile time if! how cool
      	      { //same as the last comment just with minus (reverse sign)
      		return std::make_pair(LINESTATE_LINE,(
      		    Line_2(normal.dy()*m_normalOfMainPlane.dx(),
      			normal.dz()*m_normalOfMainPlane.dx(),
			projToPlaneThroughOrigin?0:(isOnNegetive?-normal.dx():normal.dx())
			))
      		);

      	      }

      	  }
      	};
      template <typename Kernel,bool isOnNegetive>
      class PlaneProjector: public PlaneProjectorFather<Kernel>
      {
	PlaneProjectorFather<Kernel>* innerPlaneProjector;
      public:
      	virtual typename Kernel::Direction_3 point_2ToDirection_3(typename Kernel::Point_2 p)
      	{
      	return innerPlaneProjector->point_2ToDirection_3(p);
      	}
	virtual std::pair<enum LineState,typename Kernel::Line_2> projectHemisphereToPlaneAndReturnItsComplementary(typename Kernel::Direction_3 normal,bool projToPlaneThroughOrigin){
	 return innerPlaneProjector->projectHemisphereToPlaneAndReturnItsComplementary(normal,projToPlaneThroughOrigin);
	}

	PlaneProjector(typename Kernel::Direction_3 normal)
      {
	  if(normal.dz()>0)
	    {
	      innerPlaneProjector = new InnerPlaneProjector<Kernel,1,0,0,isOnNegetive>(normal);
	    }
	  else if(normal.dz()<0)
	    {
	      innerPlaneProjector = new InnerPlaneProjector<Kernel,-1,0,0,isOnNegetive>(normal);

	    }
	  else//normal.dz()==0
	    {
	      if(normal.dy()>0)
		{
		  innerPlaneProjector = new InnerPlaneProjector<Kernel,0,1,0,isOnNegetive>(normal);

		}
	      else if(normal.dy()<0)
		{
		  innerPlaneProjector = new InnerPlaneProjector<Kernel,0,-1,0,isOnNegetive>(normal);
		}
	      else //normal.dy()==0
		{
		  if(normal.dx()>0)
		    {
		      innerPlaneProjector = new InnerPlaneProjector<Kernel,0,0,1,isOnNegetive>(normal);
		    }
		  else //if(normal.dz()<0)
		    {
		      innerPlaneProjector = new InnerPlaneProjector<Kernel,0,0,-1,isOnNegetive>(normal);
		    }
		}
	    }

      }
	virtual ~PlaneProjector()
	{
	  delete innerPlaneProjector;
	}
      };
    } // end of namespace internal
  } // end of namespace Set_movable_separability_3
} // end of namespace CGAL
#endif /* INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_PLANEPROJECTOR_H_ */
