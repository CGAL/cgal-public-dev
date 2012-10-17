// Copyright (c) 1997-2012  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer
//                 Yves Brise


#ifndef CGAL_QP_DEBUG_H
#define CGAL_QP_DEBUG_H

#include <CGAL/Timer.h>
#include <CGAL/IO/Verbose_ostream.h>

// macro definitions
// =================

// debug
// -----

// make CGAL_QP_NO_DEBUG the decisive "no debug" directive
#if (    defined( CGAL_QP_NO_DEBUG)\
      || defined ( CGAL_NO_DEBUG)\
      || defined ( NDEBUG))
#  ifndef CGAL_QP_NO_DEBUG
#  define CGAL_QP_NO_DEBUG
#  endif
#endif



#if (    defined( CGAL_QP_NO_DEBUG)\
      || defined( CGAL_QP_NO_ASSERTIONS)\
      || defined( CGAL_NO_ASSERTIONS)\
      || defined( CGAL_NO_DEBUG) || defined( NDEBUG))
#  define  CGAL_qpe_debug  if ( 0)
#else
#  define  CGAL_qpe_debug  if ( 1)
#endif // qpe debug



// timing functionality
namespace CGAL {
  namespace QP_solver_debug {
    struct Timing {
      CGAL::Timer total;
      CGAL::Timer pricing;
      CGAL::Timer ratio_test;
      CGAL::Timer ratio_test_1;
      CGAL::Timer ratio_test_2;
      CGAL::Timer ratio_test_3;
      CGAL::Timer compute_factorization;
      CGAL::Timer get_basis_matrix;
      CGAL::Timer gen_1;
      CGAL::Timer gen_2;
      CGAL::Timer gen_3;
      CGAL::Timer gen_4;
      CGAL::Timer gen_5;
      
      int counter_Oin;
      int counter_Oout;
      int counter_Sin;
      int counter_Sout;
      int counter_O_O;
      int counter_O_S;
      int counter_S_O;
      int counter_S_S;
      int counter_UZ1;
      int counter_UZ2;
      int counter_UZ3;
      int counter_UZ4;
      
      int counter_Oin_fail;
      int counter_Oout_fail;
      int counter_Sin_fail;
      int counter_Sout_fail;
      int counter_O_O_fail;
      int counter_O_S_fail;
      int counter_S_O_fail;
      int counter_S_S_fail;
      int counter_UZ1_fail;
      int counter_UZ2_fail;
      int counter_UZ3_fail;
      int counter_UZ4_fail;
      
      int counter_1;
      int counter_2;
      int counter_3;
      int counter_integral_division;
      
      int matrix_size;
      int bit_size;
      
      double avg_1;
      int avg_counter_1;
      
      double density;
      
      Timing() {
        reset();
      }
      
      void reset() {
        total.reset();
        pricing.reset();
        ratio_test.reset();
        ratio_test_1.reset();
        ratio_test_2.reset();
        ratio_test_3.reset();
        compute_factorization.reset();
        get_basis_matrix.reset();
        gen_1.reset();
        gen_2.reset();
        gen_3.reset();
        gen_4.reset();
        gen_5.reset();
        
        counter_Oin = 0;
        counter_Oout = 0;
        counter_Sin = 0;
        counter_Sout = 0;
        counter_O_O = 0;
        counter_O_S = 0;
        counter_S_O = 0;
        counter_S_S = 0;
        counter_UZ1 = 0;
        counter_UZ2 = 0;
        counter_UZ3 = 0;
        counter_UZ4 = 0;
        
        counter_Oin_fail = 0;
        counter_Oout_fail = 0;
        counter_Sin_fail = 0;
        counter_Sout_fail = 0;
        counter_O_O_fail = 0;
        counter_O_S_fail = 0;
        counter_S_O_fail = 0;
        counter_S_S_fail = 0;
        counter_UZ1_fail = 0;
        counter_UZ2_fail = 0;
        counter_UZ3_fail = 0;
        counter_UZ4_fail = 0;
        
        counter_1 = 0;
        counter_2 = 0;
        counter_3 = 0;
        
        counter_integral_division = 0;
        
        matrix_size = 0;
        bit_size = 0;
        
        avg_1 = 0.0;
        avg_counter_1 = 0;
        
        density = 0.0;
      }
      
      void print(Verbose_ostream vout) {
        vout << "============\nTiming output\n============\n";
        vout << "Total: " << total.time() << "\n";
        vout << "Pricing: " << pricing.time() << "\n";
        vout << "Ratio test total: " << ratio_test.time() << "\n";
        vout << "Ratio test 1: " << ratio_test_1.time() << "\n";
        vout << "Ratio test 2: " << ratio_test_2.time() << "\n";
        vout << "Ratio test 3: " << ratio_test_3.time() << "\n";
        vout << "Gen 1: " << gen_1.time() << "\n";
        vout << "Gen 2: " << gen_2.time() << "\n";
        vout << "Gen 3: " << gen_3.time() << "\n";
        //vout << "Gen 4: " << gen_4.time() << "\n";
        //vout << "Gen 5: " << gen_5.time() << "\n";
        vout << "Factorization: " << compute_factorization.time() << "\n";
        vout << "get_basis_matrix: " << get_basis_matrix.time() << "\n";
        vout << "Counter O in: " << "Total: " << (counter_Oin+counter_Oin_fail) << ", success rate: "
             << static_cast<double>(counter_Oin)/(counter_Oin+counter_Oin_fail) << std::endl;
        vout << "Counter O out: " << "Total: " << (counter_Oout+counter_Oout_fail) << ", success rate: "
             << static_cast<double>(counter_Oout)/(counter_Oout+counter_Oout_fail) << std::endl;
        vout << "Counter S in: " << "Total: " << (counter_Sin+counter_Sin_fail) << ", success rate: "
             << static_cast<double>(counter_Sin)/(counter_Sin+counter_Sin_fail) << std::endl;
        vout << "Counter S out: " << "Total: " << (counter_Sout+counter_Sout_fail) << ", success rate: "
             << static_cast<double>(counter_Sout)/(counter_Sout+counter_Sout_fail) << std::endl;
        vout << "Counter O/O: " << "Total: " << (counter_O_O+counter_O_O_fail) << ", success rate: "
             << static_cast<double>(counter_O_O)/(counter_O_O+counter_O_O_fail) << std::endl;
        vout << "Counter O/S: " << "Total: " << (counter_O_S+counter_O_S_fail) << ", success rate: "
             << static_cast<double>(counter_O_S)/(counter_O_S+counter_O_S_fail) << std::endl;
        vout << "Counter S/O: " << "Total: " << (counter_S_O+counter_S_O_fail) << ", success rate: "
             << static_cast<double>(counter_S_O)/(counter_S_O+counter_S_O_fail) << std::endl;
        vout << "Counter S/S: " << "Total: " << (counter_S_S+counter_S_S_fail) << ", success rate: "
             << static_cast<double>(counter_S_S)/(counter_S_S+counter_S_S_fail) << std::endl;
        vout << "Counter UZ1: " << "Total: " << (counter_UZ1+counter_UZ1_fail) << ", success rate: "
             << static_cast<double>(counter_UZ1)/(counter_UZ1+counter_UZ1_fail) << std::endl;
        vout << "Counter UZ2: " << "Total: " << (counter_UZ2+counter_UZ2_fail) << ", success rate: "
             << static_cast<double>(counter_UZ2)/(counter_UZ2+counter_UZ2_fail) << std::endl;
        vout << "Counter UZ3: " << "Total: " << (counter_UZ3+counter_UZ3_fail) << ", success rate: "
             << static_cast<double>(counter_UZ3)/(counter_UZ3+counter_UZ3_fail) << std::endl;
        vout << "Counter UZ4: " << "Total: " << (counter_UZ4+counter_UZ4_fail) << ", success rate: "
             << static_cast<double>(counter_UZ4)/(counter_UZ4+counter_UZ4_fail) << std::endl;
        vout << "Counter 1: " << counter_1 << "\n";
        vout << "Counter 2: " << counter_2 << "\n";
        vout << "Counter 3: " << counter_3 << "\n";
        vout << "#Integral Divisions: " << counter_integral_division << "\n";
        vout << "Matrix Size: " << matrix_size << "\n";
        vout << "Bit Size: " << bit_size << "\n";
        //vout << "Avg 1:" << avg_1/avg_counter_1 << ", of " << avg_counter_1 << " samples\n";
        vout << "Density: " << density << "\n";
        vout << "============" << std::endl;
      }
    };
    
    static Timing timer;
  }
}


#endif // CGAL_QP_DEBUG_H


// ===== EOF ==================================================================
