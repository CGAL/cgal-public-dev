// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_SIMPLE_RS2_CALLS_H
#define CGAL_RS_SIMPLE_RS2_CALLS_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <rs_exports.h>

#ifdef CGAL_RS_OLD_INCLUDES
#define CGALRS_PTR(a)   long int a
#else
#define CGALRS_PTR(a)   void *a
#endif

namespace CGAL{
namespace RS2{

struct Simple_rs2_calls{

        static void init_solver(){
                static bool first=true;
                if(first){
                        first=false;
                        rs_init_rs();
                }else
                        rs_reset_all();
        }

        static void create_rs_upoly(mpz_t *poly,
                                    const int deg,
                                    CGALRS_PTR(ident_pol)){
                CGALRS_PTR(ident_mon);
                CGALRS_PTR(ident_coeff);
                rs_import_uppring((char*)"T");
                for(int i=0;i<=deg;++i)
                        if(mpz_sgn(poly[i])){   // don't add if == 0
                                ident_mon=rs_export_new_mon_upp_bz();
                                ident_coeff=rs_export_new_gmp();
                                rs_import_bz_gmp
                                        (ident_coeff,TO_RSPTR_IN(&(poly[i])));
                                rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
                                rs_dappend_list_mon_upp_bz(ident_pol,
                                                           ident_mon);
                        }
        }

        static int affiche_sols_eqs(mpfi_ptr *x){
                CGALRS_PTR(ident_sols_eqs);
                CGALRS_PTR(ident_node);
                CGALRS_PTR(ident_vect);
                CGALRS_PTR(ident_elt);
                int nb_elts;
                ident_sols_eqs=rs_get_default_sols_eqs();
                nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
                ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
                for(int i=0;i<nb_elts;++i){
                        ident_vect=rs_export_list_vect_ibfr_monnode
                                (ident_node);
                        CGAL_assertion_msg(rs_export_dim_vect_ibfr
                                                (ident_vect)==1,
                                           "vector dimension must be 1");
                        ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                        x[i]=(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                        ident_node=rs_export_list_vect_ibfr_nextnode
                                (ident_node);
                }
                return nb_elts;
        }

}; // struct Simple_rs2_calls

} // namespace RS2
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_RS2_CALLS_H
