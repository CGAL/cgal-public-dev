// Copyright (c) 2006, 2007, 2008, 2009, 2012 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_GEOTOP_LINE_BUILDER
#define CGAL_ALGEBRAIC_KERNEL_D_GEOTOP_LINE_BUILDER 1

#include <CGAL/config.h>

#include <boost/optional.hpp>
#include <boost/interprocess/smart_ptr/unique_ptr.hpp>

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <CGAL/Algebraic_structure_traits.h>

#include <CGAL/Algebraic_kernel_d/algebraic_curve_kernel_2_tools.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>
#include <CGAL/Algebraic_kernel_d/exceptions.h>

#include <CGAL/Algebraic_kernel_d/Curve_analysis_2_geotop_lifter.h>

#include <boost/numeric/interval.hpp>
#include <algorithm>
#include <utility>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4290)
#endif

#ifndef Bisolve_telemetry_code
# define Bisolve_telemetry_code(x) 
#endif

namespace CGAL {

namespace internal {

/*!
 * \brief Constructs Vert_line-objects for an algebraic curve.
 */
template<typename AlgebraicKernelWithAnalysis_2>
class Geotop_line_builder {

public:
  
    //! this instance's template parameter
    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

    typedef typename Algebraic_kernel_with_analysis_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1;

    //! The curve class.
    typedef typename Algebraic_kernel_with_analysis_2::Curve_analysis_2 Curve_analysis_2;

    //! Univariate polynomials
    typedef typename Algebraic_kernel_with_analysis_2::Polynomial_1 Polynomial_1;

    //! Bivariate polynomials
    typedef typename Algebraic_kernel_with_analysis_2::Polynomial_2 Polynomial_2;

    //! Type for Polynomial traits
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    //! The type for coefficients
    typedef typename Algebraic_kernel_with_analysis_2::Coefficient Coefficient;

    //! Type for x-values
    typedef typename Algebraic_kernel_with_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! The type for rational x-coordinates and for interval boundaries
    typedef typename Algebraic_kernel_with_analysis_2::Bound Bound;
    
    //! \brief \c Vert_line specification for critical x-values 
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    //! Type for multiplicities
    typedef typename Algebraic_kernel_with_analysis_2::Multiplicity_type Multiplicity_type;

protected:
    //! Typedef for the Interval type
    typedef boost::numeric::interval<Bound> Interval;
    
public:
    //! approximate coefficients 
    typedef CGAL::internal::Bitstream_coefficient_kernel_at_alpha< Algebraic_kernel_d_1 > 
      Bitstream_coefficient_kernel;           
    
    //! the base isolator used
    typedef CGAL::internal::Geotop_lifter< Bitstream_coefficient_kernel >
      Status_line_isolator;

protected:
    
    //! bi diff vanish type
    typedef typename Status_line_isolator::Bi_diff_vanish_2 Bi_diff_vanish_2;
    

public:
    
    //!\name Constructors
    //!@{

    //! Default Constructor
    Geotop_line_builder() {}

    /*!
     * \brief Constructs the builder for the \c curve object.
     *
     * Apart from the curve itself a polynomial is passed which is expected
     * to be the primitive part of the curve.
     */
    Geotop_line_builder(Algebraic_kernel_with_analysis_2* kernel,
                        const Curve_analysis_2& curve) :
      _m_kernel(kernel), _m_curve(curve)
    {}

    //!@}

    //!\name Creations
    //!@{

 private:
    
    struct Factor_info {
      
      typename Polynomial_2::const_iterator _m_ppr_end;
      bool _m_lcoeff_vanish;
      Multiplicity_type _m_mult;
      Multiplicity_type _m_mult_res_fx_fy;
      bool _m_skip_teissier;
      bool _m_no_backup;
      
    };
    
#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
    typedef std::map< Polynomial_1, Factor_info > Factor_infos;

    Factor_infos _m_factor_infos;
#endif
    
 public:
    


    /*! 
     * \brief Creates an event line at position \c alpha for the specified 
     * curve.
     *
     * Additionally, the \c id of the event line to be created has to be
     * specfied, and
     * the number of arcs that are entering from the left and leaving to the
     * right are needed. Furthermore, the flag \c root_of_resultant tells
     * whether \c alpha is a root of the resultant of the specified curve, and
     * \c root_of_content indicates whether \c alpha is a root of the content,
     * which is equivalent to the existence of a vertical line component.
     *
     * TODO some more information on implementation using Teissier and 
     * when it might fail (exceptions & handling)
     */
    Status_line_1
    create_event_line(int id, const Algebraic_real_1& alpha, 
                      int arcs_left, int arcs_right,
                      bool root_of_resultant, bool root_of_content, 
                      int mult) 
#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
        throw(CGAL::internal::Non_generic_position_exception)
#endif
    {
     
#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "Create_event_line at x = " 
                           << CGAL::to_double(alpha)
                           << " with mult = " << mult 
                           << std::flush;
#endif

      Bisolve_telemetry_code(t_cel[1].start();)
      
      const Polynomial_2 &pp = _m_curve.primitive_polynomial_2();

      // as _m_curve.primitive_polynomial_2() is a primitive polynomial, no vertical line can occur,
      // or more precisely, it is already indicated by root_of_content == true;
      // see below

      Factor_info fi;

#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
      typename Factor_infos::iterator fit = _m_factor_infos.find(alpha.polynomial());
      if (fit == _m_factor_infos.end()) {
#endif

        Bisolve_telemetry_code(t_cel_vl.start();)
          
        // TODO use PT_2 for polynomial iterators
        fi._m_ppr_end = 
          CGAL::internal::poly_end_non_vanish_leading_term(kernel(), pp, alpha);
        
        Bisolve_telemetry_code(t_cel_vl.stop();)
          
        fi._m_mult_res_fx_fy = 0;
        
        // obtain mult_res_fx_fy
#if  CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
        fi._m_mult = mult;
        Bisolve_telemetry_code(t_lcoeff.start();)
        fi._m_lcoeff_vanish = (kernel()->sign_at_1_object()(CGAL::leading_coefficient(pp),alpha) == CGAL::ZERO);
        Bisolve_telemetry_code(t_lcoeff.stop();)
        Bisolve_out("lcoeff: "<< fi._m_lcoeff_vanish <<std::endl)

        Bisolve_telemetry_code(t_mult_rxy.start();)

        if (fi._m_lcoeff_vanish) {
          fi._m_mult = 0;
        }
        if (mult > 1) /* TODO other tests to check for singularity */ {
          if (!fi._m_lcoeff_vanish) {
            fi._m_mult_res_fx_fy = _m_curve.multiplicity_as_root_of_resultant_fx_fy(alpha);
          }
        }
        Bisolve_telemetry_code(t_mult_rxy.stop();)
#endif
     
        fi._m_skip_teissier = false;
        fi._m_no_backup = false;

#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
        fit = _m_factor_infos.insert(fit, std::make_pair(alpha.polynomial(), fi));
      } else {
        fi = fit->second;
      }
#endif

      Bisolve_telemetry_code(t_ffy_isol[1].start();)

      Bitstream_coefficient_kernel bck(kernel(), alpha);

      Bisolve_telemetry_code(c_status_lines_total[1]++;)

      // TODO timer for isolators-constructors

      // TODO 2012 use smart ptr
      Status_line_isolator *isolator;

      // TODO use PT_2 for begin/end
#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER

#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
      if (!fi._m_skip_teissier) {
#endif
        Bisolve_telemetry_code(t_ffyi_tes_isol[0].start();)
        Bisolve_out("Run FastLift forever for this root of irreducible factor: " << fi._m_no_backup)
        isolator = new Status_line_isolator(pp.begin(), fi._m_ppr_end, 
                                            bck,
                                            fi._m_mult, fi._m_mult_res_fx_fy,
                                            fi._m_no_backup);
        Bisolve_telemetry_code(t_ffyi_tes_isol[0].stop();)
        Bisolve_telemetry_code(t_ffyi_tes_isol[1].start();)
        isolator->isolate();
        Bisolve_telemetry_code(t_ffyi_tes_isol[1].stop();)
#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
      } else {
        Bisolve_out("Skip FastLift for other root of irreducible factor")
      }
#endif // CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER


      if (!isolator->is_isolated()) {
        Bisolve_out("Teissier-e-fiber failed, " << std::flush);
#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
        fit->second._m_skip_teissier = true;
#endif
#endif        
#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
          delete isolator;
#endif
        if (alpha.is_rational()) {
          Bisolve_out("use isolator for exact alpha" << std::endl);
          isolator = new Status_line_isolator(bck, pp.begin(), fi._m_ppr_end, alpha.low());
          Bisolve_telemetry_code(t_ffyi_ak1_isol.start();)
          isolator->isolate();
          Bisolve_telemetry_code(t_ffyi_ak1_isol.stop();)
        } else {
          Bisolve_out("use bdv-backup" << std::endl);
          isolator = new Status_line_isolator(bi_diff_vanish(), bck, mult,
                                          // the following enables local for Teissier
                                          CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER);
          Bisolve_telemetry_code(t_ffyi_bdv_isol.start();)
          isolator->isolate();
          Bisolve_telemetry_code(t_ffyi_bdv_isol.stop();)
        }
        
          
#if  CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
      } else {
        Bisolve_telemetry_code(c_status_lines_tes[1]++;)
#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
        if (!fit->second._m_no_backup) {
          fit->second._m_no_backup = true;
          Bisolve_out("Do not stop FastLift for all other roots of irreducible factor")
        }
#endif
      }
#endif

      CGAL_assertion(isolator->is_isolated());
      
      Bisolve_telemetry_code(t_ffy_isol[1].stop();)

      Bisolve_telemetry_code(t_sllift[1].start();)

      Bisolve_telemetry_code(t_sllift_init.start();) 

     // now start to construct status line:
      int n = isolator->number_of_real_roots();
      
#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "#events: " << n << std::endl;
#endif

      std::list< int > critical;
      // copy all clusters to critical that have multiplicity > 1
      for (int j = 0; j < n; j++) {
        if (isolator->upper_bound_for_multiplicity(j) > 1) {
          critical.push_back(j);
        }
      }
      
      const int num_crit = critical.size();

#if CGAL_ACK_DEBUG_FLAG
      CGAL_ACK_DEBUG_PRINT << "#critical: " << num_crit << std::endl;
#endif
      
      typename Status_line_1::Arc_container arc_container;
      
      std::pair< int, int > minf_arcs, pinf_arcs;
      
      Bisolve_telemetry_code(t_sllift_init.stop();)

      if (num_crit == 1 && fi._m_ppr_end == pp.end()) // generic case: one critical with full degree  
      {
        Bisolve_telemetry_code(t_sllift_generic.start();)
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Generic-position case" << num_crit << std::endl;
#endif

        // critical.front() is an iterator in isolator!
        int c = critical.front();
        
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "arg(single critical): " << c << std::endl;
#endif
        
        int arcs_to_critical_left = arcs_left - n + 1;
        int arcs_to_critical_right = arcs_right - n + 1;
        
        for (int i = 0; i < n; i++) {
          if (i != c) {
            arc_container.push_back(std::make_pair(1,1));
          } else {
            arc_container.push_back(std::make_pair(arcs_to_critical_left,
                                                   arcs_to_critical_right));
          }
        }
            
        minf_arcs = pinf_arcs = std::make_pair(0,0);
        
#if 0 && !CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES
        // TODO 2012: check with MKMK; activate this, or let it done by bucketing?
        if (kernel()->is_zero_at_1_object() 
            (CGAL::leading_coefficient(pp),alpha)) {
          int d = CGAL::degree(pp,1);
          CGAL_assertion(! kernel()->is_zero_at_1_object()(CGAL::get_coefficient(pp,d-1),alpha));
          CGAL::Sign asym_sign 
            = kernel()->sign_at_1_object()
            (CGAL::get_coefficient(pp,d-1),alpha)
            * kernel()->sign_at_1_object()
            (CGAL::differentiate
             (CGAL::get_coefficient(pp,d)),alpha);
          CGAL_assertion(asym_sign!=CGAL::ZERO);
          if (asym_sign==CGAL::SMALLER) {
            minf_arcs = std::make_pair(1,0); 
            pinf_arcs = std::make_pair(0,1);
          } else {
            minf_arcs = std::make_pair(0,1); 
            pinf_arcs = std::make_pair(1,0);
          }
        }
#endif
       
        Bisolve_telemetry_code(t_sllift_generic.stop();) 
      } else {
        
        // do bucketing
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Non-Generic-position case, bucketing" << std::endl;
#endif
        
        // Remark: The following is borrowed from Michael Kerber's experimental code in Curve_analysis_2

        Bisolve_telemetry_code(t_sllift_bucket.start();)        

        // Now adjacencies
        std::vector<Bound> bucket_borders;
        
        if (n == 0) {
          bucket_borders.push_back(0);
        } else {
          bucket_borders.push_back(
                                   CGAL::internal::bound_left_of
                                   (kernel(),Algebraic_real_1(isolator->left_bound(0) -1)));
          for (int i = 1; i < n; i++) {
            while(Algebraic_real_1(isolator->right_bound(i-1))==
                  Algebraic_real_1(isolator->left_bound(i))) {
              isolator->refine_interval(i-1);
              isolator->refine_interval(i);
            }
            bucket_borders.push_back(
                                     kernel()->bound_between_1_object()
                                     (Algebraic_real_1(isolator->right_bound(i-1)),
                                      Algebraic_real_1(isolator->left_bound(i)))
                                     );
          }
          
          bucket_borders.push_back(
                                   CGAL::internal::bound_right_of
                                   (kernel(),
                                    Algebraic_real_1(isolator->right_bound(n-1) + 1)));
        }

        /*
        for (int i = 0; i <= n; i++) {
        
          std::cout << "Bucket " << i << " = " << CGAL::to_double(bucket_borders[i]) << std::endl;
          //std::cout << "Bucket " << i << " = " << (bucket_borders[i]) << std::endl;
        }
        */


        Bisolve_telemetry_code(t_sllift_bucket.stop();)        

        Bisolve_telemetry_code(t_sllift_refine.start();)

        Bound left = _m_curve.bound_value_in_interval(id);
        Bound right = _m_curve.bound_value_in_interval(id+1);
        
        typedef CGAL::Coercion_traits<Bound, Coefficient> Coercion;

        typedef typename Coercion::Type Coercion_type;

        typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
          ::template Rebind<Coercion_type,1>::Other::Type Poly_coer_1;
        
        typedef boost::numeric::interval<Coercion_type> Coercion_interval;

        for (int i = 0; i < static_cast<int>(bucket_borders.size()); i++) {
            
          Poly_coer_1 curr_pol =
            // TODO use PT_d?
            pp.evaluate(bucket_borders[i]);

          CGAL::internal::Interval_evaluate_1< Poly_coer_1, Bound >
            interval_evaluate_1;

          while(true) {

              std::pair<Bound,Bound> curr_interval_pair =
              interval_evaluate_1(curr_pol,std::make_pair(left,right));

            Coercion_interval curr_interval(curr_interval_pair.first,
                                            curr_interval_pair.second);            

            if (boost::numeric::in_zero(curr_interval)) {
              // "refine"
              Bound middle = (left+right)/2;
              CGAL::Comparison_result cmp = CGAL::compare(Algebraic_real_1(middle), alpha);
              if (cmp == CGAL::EQUAL) {
                left  = (left + middle)/2;
                right = (right + middle)/2;
              } else if (cmp == CGAL::LARGER) {
                right = middle;
              } else {
                left = middle;
              }
            } else {
              break;
            }
          }
        }

        // std::cout << "bleft: " << left << std::endl;
        // std::cout << "bright: " << right << std::endl;

        // std::cout << "#left: " << arcs_left << std::endl;
        // std::cout << "#right: " << arcs_right << std::endl;

        Bisolve_telemetry_code(t_sllift_refine.stop();)


        Bisolve_telemetry_code(t_sllift[1].stop();)

        Status_line_1 left_line =
          _m_curve.status_line_at_exact_x(Algebraic_real_1(left));
        Status_line_1 right_line =
          _m_curve.status_line_at_exact_x(Algebraic_real_1(right));
        
        Bisolve_telemetry_code(t_sllift[1].start();)

        Bisolve_telemetry_code(t_sllift_assign.start();)

        // std::cout << "leftl: " << left_line << std::endl;
        // std::cout << "rightl: " << right_line << std::endl;

        std::vector<int> left_arcs(bucket_borders.size()+1);
        std::vector< int > right_arcs(bucket_borders.size()+1);
        
        for (unsigned int i = 0; i < left_arcs.size(); i++) {
          left_arcs[i]=0;
        }
        for (unsigned int i = 0; i < right_arcs.size(); i++) {
          right_arcs[i]=0;
        }

        int curr_index = 0;
        for(int i = 0; i < arcs_left; i++) {
          while (true) {
            if (curr_index == static_cast<int>(bucket_borders.size())) {
              left_arcs[curr_index]++;
              break;
            } else if (left_line.lower_bound(i) > bucket_borders[curr_index]) {
              curr_index++;
            } else if (left_line.upper_bound(i) < bucket_borders[curr_index]) {
              left_arcs[curr_index]++;
              break;
            } else {
              left_line.refine(i);
            }
          }
        }

        curr_index = 0;
        for (int i = 0; i < arcs_right; i++) {
          while (true) {
            if( curr_index == static_cast<int>(bucket_borders.size())) {
              right_arcs[curr_index]++;
              break;
            } else if (right_line.lower_bound(i) > bucket_borders[curr_index]) {
              curr_index++;
            } else if (right_line.upper_bound(i) < bucket_borders[curr_index]) {
              right_arcs[curr_index]++;
              break;
            } else {
              right_line.refine(i);
            }
          }
        }

        for (int i = 0; i < n; i++) {
          arc_container.push_back(std::make_pair(left_arcs[i+1],right_arcs[i+1]));
        }

        minf_arcs = std::make_pair(left_arcs[0],right_arcs[0]);
        pinf_arcs = std::make_pair(left_arcs[n+1],right_arcs[n+1]);

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

        Bisolve_telemetry_code(t_sllift_assign.stop();)

      }
      
      Bisolve_telemetry_code(t_sllift_finish.start();)

      Status_line_1 vl(alpha, id, _m_curve, arcs_left, arcs_right, arc_container);
      
      vl._set_number_of_branches_approaching_infinity(minf_arcs, pinf_arcs);
      
      vl.set_isolator(*isolator);
      
      // set vertical line
      if (root_of_content) {
        vl._set_v_line();
      }

      Bisolve_telemetry_code(t_sllift_finish.stop();)

      Bisolve_telemetry_code(t_sllift[1].stop();)

      Bisolve_telemetry_code(t_cel[1].stop();)

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done CESL " 
                             << vl
                             << std::endl;
#endif
      return vl;
      
    }
    
    /*! 
     * \brief Creates a non-event line at position \c ar for the specified 
     * curve.
     *
     * Additionally, the \c id of the event line to be created has to be
     * specfied.
     *
     * TODO some more information on implementation using Teissier and 
     * when it might fail (exceptions & handling)
     */
    Status_line_1
    create_non_event_line(int id, const Algebraic_real_1& ar) {
      
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Create_non_event_line at x = " 
                             << CGAL::to_double(ar)
                             << " (" << ar << ")"
                             << std::flush;
#endif
        
      Bisolve_telemetry_code(t_cel[2].start();)

      const Polynomial_2& pp = _m_curve.primitive_polynomial_2(); // must have full degree as called for ar != "event"
      
      //std::cout << "pp: " << pp << std::endl;
      //std::cout << "ar: " << ar << std::endl;
      
      Bisolve_telemetry_code(t_ffy_isol[2].start();)

      Bitstream_coefficient_kernel bck(kernel(), ar);
     
      Bisolve_telemetry_code(c_status_lines_total[2]++;)
      
      // TODO timer for isolators-constructors
      
        // TODO 2012 use unique_ptr? optional? avoid copy and assigment
        // Status_line_isolator isolator(_m_curve.polynomial_2().begin(), _m_curve.polynomial_2().end(), bck); // dummy construction

        Status_line_isolator *isolator = 0;

      // TODO use PT_2 for iterators to pp
#if  0 && CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
      isolator = new Status_line_isolator(bck, pp.begin(), pp.end(), 0, 0);
      Bisolve_telemetry_code(t_ffyi_tes_isol[2].start();)
      isolator->isolate();
      Bisolve_telemetry_code(t_ffyi_tes_isol[2].stop();)
      if (!isolator->is_isolated()) {
        std::cout << "Teissier-ne-fiber failed, " << std::flush;
#endif
#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
          delete isolator;
#endif
        if (ar.is_rational()) {
          Bisolve_out("use isolator for exact alpha backup" << std::endl);
          isolator = new Status_line_isolator(bck, pp.begin(), pp.end(), ar.low());
        } else {
          Bisolve_out("use bitstream backup" << std::endl);
          isolator = new Status_line_isolator(pp, bck);

        }
        Bisolve_telemetry_code(t_ffyi_bits_isol.start();)
        isolator->isolate();
        Bisolve_telemetry_code(t_ffyi_bits_isol.stop();)

#if 0 && CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
      } else {
        Bisolve_telemetry_code(c_status_lines_tes[2]++;)
      }
#endif      

      CGAL_assertion(isolator->is_isolated());

      Bisolve_telemetry_code(t_ffy_isol[2].stop();)

      Bisolve_telemetry_code(t_sllift[2].start();)

      // now start to construct status line:
      size_t root_number = isolator->number_of_real_roots();
      
      //std::cout << "rn: " << root_number << std::endl;
      
      Status_line_1 status_line(ar, id, _m_curve, root_number);
      
      status_line.set_isolator(*isolator);
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done CNESL " 
                             << status_line
                             << std::endl;
#endif
      Bisolve_telemetry_code(t_sllift[2].stop();)

      Bisolve_telemetry_code(t_cel[2].stop();)

      return status_line;
    }
    
    //!@}

protected:

    //! the univariate kernel
    Algebraic_kernel_with_analysis_2* kernel() const {
        return this->_m_kernel;
    }

public:

    //! returns bdv instance
    const Bi_diff_vanish_2& bi_diff_vanish() const {

      if (!_m_bdv) {
        _m_bdv = Bi_diff_vanish_2(_m_curve.primitive_polynomial_2());
      }
      
      return *_m_bdv;
    }
    

protected:
    
    //! the univariate kernel
    Algebraic_kernel_with_analysis_2* _m_kernel;

    //! The curve whose Status_line_1s are built.
    Curve_analysis_2 _m_curve;

    //! map for bi diff vanishs
    mutable boost::optional< Bi_diff_vanish_2 > _m_bdv;

}; //class Geotop_line_builder

} // namespace internal

} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_ALGEBRAIC_KERNEL_D_GEOTOP_LINE_BUILDER
// EOF
