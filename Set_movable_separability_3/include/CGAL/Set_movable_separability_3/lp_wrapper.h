/*
 * lp_wrapper.h
 *
 *  Created on: Apr 8, 2017
 *      Author: shahar
 */

#ifndef SET_MOVABLE_SEPARABILITY_3_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_LP_WRAPPER_H_
#define SET_MOVABLE_SEPARABILITY_3_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_LP_WRAPPER_H_
/*
 * This function should be redundant.
 * In our solution we use LP multiple times (in 2D).
 * In one of this times, we know that the intersection is empty and
 * we need this the LP to return 3 constrains (half-planes) that their intersection is empty.
 * Such 3 constrains must exist according to "Helly's theorem" - full proof is in our paper:
 * "On the Separation of a Polyhedron from Its Single-Part Mold" (version of 2017).
 *
 * Sadly the LP package in CGAL "Linear and Quadratic Programming Solver" do not offer such interface.
 * According to Bernd Gärtner this can be deduced from the package "infeasibility_certificate" thanks
 * to the current implementation.
 *
 * This function convert the "infeasibility_certificate" to this 3 half-planes.
 * If it fails (it might fail if "Linear and Quadratic Programming Solver" implementation will change),
 * An Exception will be thrown.
 */

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
namespace CGAL {
namespace Set_movable_separability_3 {
namespace internal {
//

    enum HalfplaneInteractionState
    {
      HalfplaneInteractionState_LinesIntersect,
      HalfplaneInteractionState_LinesParallel_FisrtContainsSecond, //equal is in one of this
      HalfplaneInteractionState_LinesParallel_SecondContainsFirst, //equal is in one of this
      HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty,
      HalfplaneInteractionState_LinesParallel_UnionIsAll //not possible under our conditions,
    };
    /*
    * a1 b1 c1  	 	=> 	a1*x+b1*y <= c1
    * a2 b2 c2  		=> 	a2*x+b2*y <= c2
    */
    template <typename FT>
    enum HalfplaneInteractionState halfPlaneInteraction(FT *h1,FT* h2)
    {
	if(h1[0]*h2[1]!=h2[0]*h1[1])
	{
	   return HalfplaneInteractionState_LinesIntersect;
	}
	if(h1[0]!=0)
	{
	    bool h1xPositive=(h1[0]>0);
	    bool h2xPositive=(h2[0]>0);
	    if((h1xPositive&&(!h2xPositive))||(h2xPositive&&(!h1xPositive)))
	    {
	      return HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty; //HalfplaneInteractionState_LinesParallel_UnionIsAll might be as well, but In our case, the input to this function can't be in this state
	    }
	}
	else
	{
	    bool h1yPositive=(h1[1]>0);
	    bool h2yPositive=(h2[1]>0);
	    if((h1yPositive&&(!h2yPositive))||(h2yPositive&&(!h1yPositive)))
	      {
		return HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty; //HalfplaneInteractionState_LinesParallel_UnionIsAll might be as well, but In our case, the input to this function can't be in this state
	      }
	}
	if(h1[2]>h2[0])
	   return HalfplaneInteractionState_LinesParallel_SecondContainsFirst;
	else
	   return HalfplaneInteractionState_LinesParallel_FisrtContainsSecond;

    }

   /*this function solves finds a witnesses of an infeasiblity of the input lp.
   * (pre: this lp is infeasible) and return 3 (or 2) half planes proving its infeasiblity
   * input: halfplanes -  array of halfplanes in format:
   * a1 b1 c1  	 	=> 	a1*x+b1*y <= c1
   * a2 b2 c2  		=> 	a2*x+b2*y <= c2
   * a3 b3 c3   	=> 	a3*x+b3*y <= c3
   * ....
   * an bn cn   	=> 	an*x+bn*y <= cn
   *
   * and pointer to an empty array of size 3 uints at least
   *
   * output: the number of half planes inserted to "outArray" (2 or 3)
   * notice that is guaranteed that none of this half-planes are redundant! if you get 3 half-planes, you need all three of them to prove
   * the intersection is empty.
   * This do not means that if you get 3 half-planes there can't be two half-planes that could use as witnesses.
   */
   template <typename FT>
   uint8_t find3Witnesses(FT halfplanes[][3],unsigned int n, unsigned int * outArray)
   {
     typedef CGAL::Quadratic_program<FT> Program;
     typedef CGAL::Quadratic_program_solution<FT> Solution;
     Program lp(CGAL::SMALLER, false,0,false,0);
#define X_LP_INDEX 0
#define Y_LP_INDEX 1
     for(int i=0;i<n;++i)
     {
	 lp.set_a(X_LP_INDEX, i, halfplanes[i][X_LP_INDEX]);
	 lp.set_a(Y_LP_INDEX, i, halfplanes[i][Y_LP_INDEX]);
	 lp.set_b(i++, halfplanes[i][2]);
     }
     Solution s;
     s= CGAL::solve_linear_program(lp, FT());

     if(!s.is_infeasible())
     {
	 throw std::invalid_argument("program is feasibile");
     }
     uint32_t i=0;
     uint8_t count=0;
     for(typename Solution::Infeasibility_certificate_iterator it= s.infeasibility_certificate_begin(); it!= s.infeasibility_certificate_end();it++)
     {
	 if(*it!=0)
	 {
	     if(count==3)
	     {
		 throw std::invalid_argument("The Infeasibility_certificate_iterator returned more then 3 half planes! BUG!");
	     }
	    outArray[count++]=i;
	 }
	 i++;
     }
     if(count==3)
     {
	 for(int i=0;i<2;++i)
	 {
	     for(int j=i+1;j<3;++j)
	     {
		 switch (halfPlaneInteraction(halfplanes[i],halfplanes[j]))
		 {

		   case HalfplaneInteractionState_LinesParallel_FisrtContainsSecond:
		     //in this case we want to delete the first
		     outArray[i]=outArray[2];
		     return 2;
		   case HalfplaneInteractionState_LinesParallel_SecondContainsFirst:
		     //in this case we want to delete the second
		     outArray[j]=outArray[2];
		     return 2;
		   case HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty:
		     //in this case we want to delete the third halfplane
		     outArray[(0+1+2)-i-j]=outArray[2];
		     return 2;
		 }
	     }

	 }
	 return 3;
     }
     else if(count<2)
     {
	 throw std::invalid_argument("The Infeasibility_certificate_iterator returned less then 2 half planes! BUG!");
     }

     return 2;
   }
} // end of namespace internal
} // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif /* SET_MOVABLE_SEPARABILITY_3_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_LP_WRAPPER_H_ */
