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
    template <typename Kernel>
    enum HalfplaneInteractionState halfPlaneInteraction(typename Kernel::Line_2 &h1,typename Kernel::Line_2 &h2)
    {
	if(h1.a()*h2.b()!=h2.a()*h1.b())
	{
	   return HalfplaneInteractionState_LinesIntersect;
	}
	if(h1.a()!=0)
	{
	    bool h1xPositive=(h1.a()>0);
	    bool h2xPositive=(h2.a()>0);
	    if((h1xPositive&&(!h2xPositive))||(h2xPositive&&(!h1xPositive)))
	    {
	      return HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty; //HalfplaneInteractionState_LinesParallel_UnionIsAll might be as well, but In our case, the input to this function can't be in this state
	    }
	}
	else
	{
	    bool h1yPositive=(h1.b()>0);
	    bool h2yPositive=(h2.b()>0);
	    if((h1yPositive&&(!h2yPositive))||(h2yPositive&&(!h1yPositive)))
	      {
		return HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty; //HalfplaneInteractionState_LinesParallel_UnionIsAll might be as well, but In our case, the input to this function can't be in this state
	      }
	}
	if(h1.c()>h2.c())
	   return HalfplaneInteractionState_LinesParallel_SecondContainsFirst;
	else
	   return HalfplaneInteractionState_LinesParallel_FisrtContainsSecond;

    }

   /*this function solves finds a witnesses of an infeasiblity of the input lp.
   * (pre: this lp is infeasible) and return 3 (or 2) half planes proving its infeasiblity
   * input: halfplanes -  array of halfplanes in format:
   * line_2 , name as int. if the line is (a,b,c) it represents ax+by+c>=0
   *
   * and pointer to an empty array of size 3 uints at least
   *
   * output: the number of half planes inserted to "outArray" (2 or 3)
   * notice that is guaranteed that none of this half-planes are redundant! if you get 3 half-planes, you need all three of them to prove
   * the intersection is empty.
   * This do not means that if you get 3 half-planes there can't be two half-planes that could use as witnesses.
   */
   template <typename Kernel>
   uint8_t find3Witnesses(const std::vector<std::pair<typename Kernel::Line_2,unsigned int>> halfplanes, unsigned int * outArray)
   {
     typedef typename Kernel::Line_2 Line_2;
     typedef typename Kernel::FT FT;
     typedef CGAL::Quadratic_program<FT> Program;
     typedef CGAL::Quadratic_program_solution<FT> Solution;
     Program lp(CGAL::LARGER, false,0,false,0);
#define X_LP_INDEX 0
#define Y_LP_INDEX 1
     for(int i=0;i<halfplanes.size();++i)
     {
	 //Line_2 format is (a,b,c) => ax+by+c>=0, the format here is ax+by>=c
	 lp.set_a(X_LP_INDEX, i, halfplanes[i].first.a());
	 lp.set_a(Y_LP_INDEX, i, halfplanes[i].first.b());
	 lp.set_b(i++, -halfplanes[i].c());
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
		 switch (halfPlaneInteraction<Kernel>(halfplanes[outArray[i]].first,halfplanes[outArray[j]].first))
		 {

		   case HalfplaneInteractionState_LinesParallel_FisrtContainsSecond:
		     //in this case we want to delete the first
		     outArray[i]=outArray[2];
		     outArray[0]=halfplanes[outArray[0]].second;
		     outArray[1]=halfplanes[outArray[1]].second;
		     return 2;
		   case HalfplaneInteractionState_LinesParallel_SecondContainsFirst:
		     //in this case we want to delete the second
		     outArray[j]=outArray[2];
		     outArray[0]=halfplanes[outArray[0]].second;
		     outArray[1]=halfplanes[outArray[1]].second;
		     return 2;
		   case HalfplaneInteractionState_LinesParallel_Intesection_Is_Empty:
		     //in this case we want to delete the third halfplane
		     outArray[(0+1+2)-i-j]=outArray[2];
		     outArray[0]=halfplanes[outArray[0]].second;
		     outArray[1]=halfplanes[outArray[1]].second;
		     return 2;
		 }
	     }

	 }
	 outArray[0]=halfplanes[outArray[0]].second;
	 outArray[1]=halfplanes[outArray[1]].second;
	 outArray[2]=halfplanes[outArray[2]].second;
	 return 3;
     }
     else if(count<2)
     {
	 throw std::invalid_argument("The Infeasibility_certificate_iterator returned less then 2 half planes! BUG!");
     }
     outArray[0]=halfplanes[outArray[0]].second;
     outArray[1]=halfplanes[outArray[1]].second;
     return 2;
   }
} // end of namespace internal
} // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif /* SET_MOVABLE_SEPARABILITY_3_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_3_LP_WRAPPER_H_ */
