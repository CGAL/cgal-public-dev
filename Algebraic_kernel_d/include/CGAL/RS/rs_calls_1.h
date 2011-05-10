// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_RS_CALLS_1_H
#define CGAL_RS_RS_CALLS_1_H

#include <stdlib.h>
#include <CGAL/RS/basic_1.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <rs_exports.h>

namespace CGAL{

#define CGALRS_CSTR(S)  ((char*)(S))
  
#ifdef CGAL_RS_OLD_INCLUDES
#define CGALRS_PTR(a)   long int a
#else
#define CGALRS_PTR(a)   void *a
#endif
  
// initialize RS solver
inline void init_solver(){
        static bool first=true;
        if(first){
                first=false;
                //rs_version(stderr);
                rs_init_rs();
        }else{
                rs_reset_all();
        }
}

// reset RS memory
inline void reset_solver(){
        rs_reset_all();
}

inline int affiche_sols_eqs(mpfi_ptr *x){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts;
        ident_sols_eqs=rs_get_default_sols_eqs();
        nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
        ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
        //x=(mpfi_ptr*)malloc(nb_elts*sizeof(mpfi_ptr));
        for(int i=0;i<nb_elts;++i){
                ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
                CGAL_assertion_msg(rs_export_dim_vect_ibfr(ident_vect)==1,
                                "the dimension of vector must be 1");
                ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                x[i]=(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt);
                ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
        }
        return nb_elts;
}

inline void affiche_sols_constr(int nr,mpfi_ptr p){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts,nb;
        ident_sols_eqs=rs_get_default_sols_ineqs();
        nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
        ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
        for(int i=0;i<nb_elts;++i){
                ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
                nb=rs_export_dim_vect_ibfr(ident_vect);
                CGAL_assertion_msg((nb==1),
                                "the vector must contain one element");
                ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                if(i==nr){
                        mpfi_set(p,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
                        //break;
                }
                ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
        }
}

// nroot is the number of root of the algebraic number
inline Sign affiche_signs_constr(int nroot){
        CGALRS_PTR(ident_sols_eqs);
        CGALRS_PTR(ident_node);
        CGALRS_PTR(ident_vect);
        CGALRS_PTR(ident_elt);
        int nb_elts,nb;
        mpfi_t tmp;
        mpfi_init(tmp);
        ident_sols_eqs=rs_get_default_sols_ineqs();
        nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
        ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
        for(int i=1;i<nb_elts+1;++i){
                ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
                nb=rs_export_dim_vect_ibfr(ident_vect);
                CGAL_assertion_msg((nb==1),
                                "the vector must contain one element");
                ident_elt=rs_export_elt_vect_ibfr(ident_vect,0);
                if(i==nroot+1){
                        mpfi_set(tmp,(mpfi_ptr)rs_export_ibfr_mpfi(ident_elt));
                        break;
                }
                ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
        }
        //std::cerr << "\nreturned value: ";
        //mpfi_out_str(stderr,10,0,tmp);
        //std::cerr << std::endl;
        // mpfi_is_zero(tmp) doesn't work. The reason is that MPFR_SIGN in
        // the mpfi code returns 1 when applied to the left and right zeros.
        // This is not surprising because zero is signed in IEEE 754, and MPFR
        // adopts it. Nevertheless, mpfr_sgn returns 0, but mpfi doesn't use
        // it to implement mpfi_is_zero.
        // Here is the difference (from MPFR source code):
        //  define mpfr_sgn(_x)      (mpfr_zero_p(_x) ? 0 : MPFR_SIGN(_x))
        //
        if(mpfr_zero_p(&(tmp->right))&&mpfr_zero_p(&(tmp->left)))
                return ZERO;
        // the same holds for mpfi_is_pos and mpfi_is_neg
        if((mpfr_sgn(&(tmp->left))>=0)&&(mpfr_sgn(&(tmp->right)))>0)
                return POSITIVE;
        if((mpfr_sgn(&(tmp->left))<0)&&(mpfr_sgn(&(tmp->right))<=0))
                return NEGATIVE;
        // if we arrive here, it is because the signs of the endpoints are -
        // and +, and (I think) RS guarantees that this never happens
        CGAL_assertion_msg(false,"error in sign calculation");
        return ZERO;
}

inline void create_rs_upoly(mpz_t *poly,const int deg,CGALRS_PTR(ident_pol)){
        CGALRS_PTR(ident_mon);
        CGALRS_PTR(ident_coeff);
        rs_import_uppring(CGALRS_CSTR("T"));
        for(int i=0;i<=deg;++i)
                if(mpz_sgn(poly[i])){   // don't add if == 0
		  ident_mon=rs_export_new_mon_upp_bz();
		  ident_coeff=rs_export_new_gmp();
		  rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(&(poly[i])));
		  rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
		  rs_dappend_list_mon_upp_bz(ident_pol,ident_mon);
                }
}

inline void create_rs_uconstr(mpz_t **list_constr,
                               const int *list_degs,
                               CGALRS_PTR(ident_list)){
        CGALRS_PTR(ident_poly);
        ident_poly=rs_export_new_list_mon_upp_bz();
        create_rs_upoly(*list_constr,*list_degs,ident_poly);
        rs_dappend_list_sup_bz(ident_list,ident_poly);
}
//************************************************Yacine code***************
 
// create an Rs univariate polynomial from MP_INT* instead of mpz_t
 
 void create_rs_upoly2(MP_INT * poly, const int deg,CGALRS_PTR(ident_pol))
{
  int i;
  // CGALRS_PTR(ident_ring);
  rs_import_uppring(CGALRS_CSTR("T"));
  CGALRS_PTR(ident_mon);
  CGALRS_PTR(ident_coeff);
  /* demande à RS l'adresse de son polynome univarie global */
  for (i=0;i<=deg;i++) {
    ident_mon=rs_export_new_mon_upp_bz();
    /* demande a RS de creer un entier long */
    ident_coeff=rs_export_new_gmp();
    /* affecte l'entier gmp de RS avec l'entier gmp local */
    rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(&(poly[i])));
    /* affecte le monôme de RS avec l'entier local et le degré qui va bien */
    rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
    /* ajoute le nouveau monome au polynôme univarié par défaut de RS */
    rs_dappend_list_mon_upp_bz(ident_pol,ident_mon);
  }
}

 // create an RS bivariate polynomial from MP_INT**
 
 int create_rs_bipoly(MP_INT ** poly,int deg1,int deg2,CGALRS_PTR(ident_pol))
{
  int i,j;
  //CGALRS_PTR(ident_ring);
  CGALRS_PTR(ident_mon);
  CGALRS_PTR(ident_coeff);
  CGALRS_PTR(ident_pp);
  /* demande à RS l'adresse de son polynome univarie global */
  for (i=0;i<=deg1;i++) {
  
    for (j=0;j<=deg2;j++) {
      
      /* demande a RS de creer un nouveau monôme */
      if (mpz_cmp_ui(&(poly[i][j]),0)) {
	
	ident_mon=rs_export_new_mon_pp_bz();
	
	/* demande a RS de creer un entier long */
	ident_coeff=rs_export_new_gmp();
	
	/* affecte l'entier gmp de RS avec l'entier gmp local */
	rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(&(poly[i][j])));
	
	/* affecte le monôme de RS avec l'entier local et le degré qui va bien */
	ident_pp=rs_export_new_pp();
	
	rs_set_pp_deg(ident_pp,0,i);
	rs_set_pp_deg(ident_pp,1,j);
	rs_dset_mon_pp_bz(ident_mon,ident_coeff,ident_pp);
	/* ajoute le nouveau monome au polynôme univarié par défaut de RS */
	rs_dappend_list_mon_pp_bz(ident_pol,ident_mon);
	 }
    }
  }
}


 //create a bivariate system in RS using two MP_INT**
int create_rs_bisys(MP_INT **p1,int deg11, int deg12,MP_INT **p2,int deg21, int deg22)
{
  rs_init_ppring_nb(2,CGALRS_CSTR("DRL"));
  rs_import_ppring_var(0,CGALRS_CSTR("y"));
  rs_import_ppring_var(1,CGALRS_CSTR("x"));

  CGALRS_PTR(ident_sys);
  CGALRS_PTR(ident_poly);
  ident_sys=rs_get_default_grob();
  
  ident_poly=rs_export_new_list_mon_pp_bz();
  
  create_rs_bipoly(p1,deg11,deg12,ident_poly);
  
  rs_dappend_list_smp_bz(ident_sys,ident_poly);
  
  ident_poly=rs_export_new_list_mon_pp_bz();
    
  create_rs_bipoly(p2,deg21,deg22,ident_poly);
  
  rs_dappend_list_smp_bz(ident_sys,ident_poly);
  
}

// the following five functions are used to display the result of a RUR computation on our bivariate system
 
 void affiche_bz(CGALRS_PTR(bz_ident))
{
  MP_INT * z=(MP_INT *)(bz_ident);
  if (mpz_cmp_si(z,0)>0) fprintf(stderr,"+");
  mpz_out_str(stderr,10,z);
}

 void affiche_bimon_bz(CGALRS_PTR(mon_ident))
 {
   CGALRS_PTR(pp_ident);
   affiche_bz(rs_export_mon_pp_bz_coeff(mon_ident));
   pp_ident=rs_export_mon_pp_bz_pp(mon_ident);
   fprintf(stderr,"*X^%d*Y^%d",rs_export_pp_deg(pp_ident,0),rs_export_pp_deg(pp_ident,1));
}

 void affiche_bipol_bz(CGALRS_PTR(pol_ident))
{
  int i=0;
  CGALRS_PTR(f_node)=rs_export_list_mon_pp_bz_firstnode(pol_ident);
  int nb_mons=rs_export_list_mon_pp_bz_nb(pol_ident);
  for (i=0;i<nb_mons;i++) {
    affiche_bimon_bz(rs_export_list_mon_pp_bz_monnode(f_node));
    f_node=rs_export_list_mon_pp_bz_nextnode(f_node);
  }
}

 void affiche_list_bipol_bz(CGALRS_PTR(list_pol_ident))
{
  int i=0;
  CGALRS_PTR(mynode)=rs_export_list_smp_bz_firstnode(list_pol_ident);
  int nb=rs_export_list_smp_bz_nb(list_pol_ident);
  for (i=0;i<nb;i++) {
    affiche_bipol_bz(rs_export_list_smp_bz_monnode(mynode));
    printf("\n");
    mynode=rs_export_list_smp_bz_nextnode(mynode);
  }
}


void affiche_bivariate_decomp_raw()
{
  fprintf(stderr,"\n Compressed RUR : ");
  affiche_list_bipol_bz(rs_get_default_grob());
}

// return the coefficient of the polynomial (pol_ident) as vector of vector
 
 std::vector< std::vector<MP_INT> > polynomial_coeff(CGALRS_PTR(pol_ident))
   {
     std::vector< std::vector<MP_INT> > vec;
     CGALRS_PTR(f_node);
     CGALRS_PTR(mon_ident);
     CGALRS_PTR(pp_ident);
     
     int nb_mons=rs_export_list_mon_pp_bz_nb(pol_ident);
     
     f_node = rs_export_list_mon_pp_bz_firstnode(pol_ident);
     for(int i=0;i<nb_mons;i++)
       {
	 mon_ident=rs_export_list_mon_pp_bz_monnode(f_node);
	 pp_ident=rs_export_mon_pp_bz_pp(mon_ident);
	 
	 int deg_y = rs_export_pp_deg(pp_ident,1);
	 
	 MP_INT* coeff = (MP_INT*)rs_export_mon_pp_bz_coeff(mon_ident);
	 if (deg_y < vec.size()){
	   vec[deg_y].push_back(*coeff) ;
	 }
	 else{
	   vec.resize(deg_y+1);
	   vec[deg_y].push_back(*coeff) ;
	 }
	 f_node=rs_export_list_mon_pp_bz_nextnode(f_node);
	 
       }
     return vec;
   }

 
 
//this function store the resulting RURs from RS resolution in a suitable data structure.

 std::vector< std::vector< std::vector<MP_INT> > > Rurs_sys_list()
   {
     
     std::vector< std::vector< std::vector<MP_INT> > > vect;
     int k;
     
     CGALRS_PTR(list_pol_ident) = rs_get_default_grob();
     int nb_pol=rs_export_list_smp_bz_nb(list_pol_ident);
     CGALRS_PTR(pol_ident)= rs_export_list_smp_bz_firstnode(list_pol_ident);
     
     for(k=0;k<nb_pol;k++)
       {
	 vect.push_back(polynomial_coeff(rs_export_list_smp_bz_monnode(pol_ident))); 
	 pol_ident=rs_export_list_smp_bz_nextnode(pol_ident);
       }
     return vect ;
   }
 

} // namespace CGAL

#endif  // CGAL_RS_RS_CALLS_1_H
