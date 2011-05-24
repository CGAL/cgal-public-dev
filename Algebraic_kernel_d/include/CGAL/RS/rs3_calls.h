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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/rs3_calls.h $
// $Id: rs3_calls.h 63060 2011-04-20 12:13:38Z penarand $
//
// 

#ifndef CGAL_RS_RS3_CALLS_H
#define CGAL_RS_RS3_CALLS_H

#include <rs_exports.h>
#include <CGAL/RS/basic_1.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfi.h>


namespace CGAL{

  namespace RS3{

#define CGALRS_CSTR(S)  ((char*)(S))
  
#ifdef CGAL_RS_OLD_INCLUDES
#define CGALRS_PTR(a)   long int a
#else
#define CGALRS_PTR(a)   void *a
#endif
  
    typedef Polynomial<Gmpz> Polynomial_1;
    typedef Polynomial<Polynomial_1> Polynomial_2;
    typedef Rur_2<Polynomial_1> rur_2;
    
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
    
    
    // create univariate Rs polynomial from a CGAL one
    
    void create_rs_upoly(Polynomial_1 poly,CGALRS_PTR(ident_pol))
      {
	int i;
	int deg = poly.degree();
	//CGALRS_PTR(ident_ring);
	CGALRS_PTR(ident_mon);
	CGALRS_PTR(ident_coeff);
	/* demande à RS l'adresse de son polynome univarie global */
	for (i=0;i<deg+1;i++) {
    
	  ident_mon=rs_export_new_mon_upp_bz();
	  /* demande a RS de creer un entier long */
	  ident_coeff=rs_export_new_gmp();
	  /* affecte l'entier gmp de RS avec l'entier gmp local */
	  rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(poly[i].mpz()));
	  /* affecte le monôme de RS avec l'entier local et le degré qui va bien */
	  if (mpz_cmp_ui(poly[i].mpz(),0)) {
	  rs_dset_mon_upp_bz(ident_mon,ident_coeff,i);
	  /* ajoute le nouveau monome au polynôme univarié par défaut de RS */
	  rs_dappend_list_mon_upp_bz(ident_pol,ident_mon);
	  }
	}
      }

    // construct a rur as list of univariate polynomials in order to isolate it
    
    void create_rs_rur(rur_2 list_constr)
      {
	rs_import_uppring((char*)"t");
	int i;
	CGALRS_PTR(ident_list);
	CGALRS_PTR(ident_poly);
	// create a rur for a bivariate systeme
	ident_list = rs_get_default_rur_num_vars();
        // add the main polynomial of the rur
	create_rs_upoly(list_constr.get_f(), rs_get_default_rur_ext());
	// add the denominator of the rur
	create_rs_upoly(list_constr.get_g(), rs_get_default_rur_den());
	// add the polynomials of the two coordinates into a list
	ident_poly=rs_export_new_list_mon_upp_bz();
	create_rs_upoly(list_constr.get_g1(),ident_poly);
	rs_dappend_list_sup_bz(ident_list,ident_poly);
	ident_poly=rs_export_new_list_mon_upp_bz();
	create_rs_upoly(list_constr.get_g2(),ident_poly);
	rs_dappend_list_sup_bz(ident_list,ident_poly);
	
      }
    
    


 

 // create an RS bivariate polynomial from CGAL polynomial

   int create_rs_bipoly(Polynomial_2 poly,CGALRS_PTR(ident_pol))
   {
     int i,j;
     //CGALRS_PTR(ident_ring);
     CGALRS_PTR(ident_mon);
     CGALRS_PTR(ident_coeff);
     CGALRS_PTR(ident_pp);
     /* demande à RS l'adresse de son polynome univarie global */
     int deg1 = poly.degree();
    
     for (i=0;i<=deg1;i++) {
       
       int deg2 = poly[i].degree();
	 
	 for (j=0;j<=deg2;j++) {
	   
	   /* demande a RS de creer un nouveau monôme */
	   if (mpz_cmp_ui(poly[i][j].mpz(),0)) {
	     
	     ident_mon=rs_export_new_mon_pp_bz();
	     
	     /* demande a RS de creer un entier long */
	     ident_coeff=rs_export_new_gmp();
	     
	     /* affecte l'entier gmp de RS avec l'entier gmp local */
	     rs_import_bz_gmp(ident_coeff,TO_RSPTR_IN(poly[i][j].mpz()));
	     
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


 //create a bivariate system in RS using two CGAL Polynomials
   
  int create_rs_bisys(Polynomial_2 p1,Polynomial_2 p2)
  {
    rs_init_ppring_nb(2,CGALRS_CSTR("LEX"));
    rs_import_ppring_var(0,CGALRS_CSTR("y"));
    rs_import_ppring_var(1,CGALRS_CSTR("x"));
  
    
  CGALRS_PTR(ident_sys);
  CGALRS_PTR(ident_poly);
  
  ident_sys=rs_get_default_grob();
  
  ident_poly=rs_export_new_list_mon_pp_bz();
  
  create_rs_bipoly(p1,ident_poly);
  
  rs_dappend_list_smp_bz(ident_sys,ident_poly);
  
  ident_poly=rs_export_new_list_mon_pp_bz();
    
  create_rs_bipoly(p2,ident_poly);
  
  rs_dappend_list_smp_bz(ident_sys,ident_poly);
  
}

// the following five functions are used to display the result of the RUR computation on our bivariate system
 
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
    printf("  \n");
    printf("  \n");
    mynode=rs_export_list_smp_bz_nextnode(mynode);
  }
}


void affiche_bivariate_decomp_raw()
{
  fprintf(stderr,"\n Compressed RUR : ");
  affiche_list_bipol_bz(rs_get_default_grob());
}

// return the coefficients of the RS polynomial (pol_ident) writeen in a monomiale form as a vector of vector of Gmpz
 
std::vector< std::vector<Gmpz> > polynomial_coeff(CGALRS_PTR(pol_ident))
   {
     std::vector< std::vector<Gmpz> > vec;
     CGALRS_PTR(f_node);
     CGALRS_PTR(mon_ident);
     CGALRS_PTR(pp_ident);
     
     int nb_mons=rs_export_list_mon_pp_bz_nb(pol_ident);
     
     f_node = rs_export_list_mon_pp_bz_firstnode(pol_ident);
     for(int i=0;i<nb_mons;i++)
       {
	 mon_ident=rs_export_list_mon_pp_bz_monnode(f_node);
	 pp_ident=rs_export_mon_pp_bz_pp(mon_ident);
	 
	 int deg_y = rs_export_pp_deg(pp_ident,0);
	 int deg_x = rs_export_pp_deg(pp_ident,1);	 

	 MP_INT* coeff = (MP_INT*)rs_export_mon_pp_bz_coeff(mon_ident);
	 CGAL::Gmpz g=CGAL::Gmpz(coeff);
	 if (deg_y < vec.size()){
	   if (deg_x <= vec[deg_y].size()) {
	     vec[deg_y].push_back(g) ;
	   }
	   else{
	     vec[deg_y].resize(deg_x);
	     vec[deg_y].push_back(g) ;
	   }
	 }
	 else{
	   vec.resize(deg_y+1);
	   if (deg_x <= vec[deg_y].size()) {
	     vec[deg_y].push_back(g) ;
	   }
	   else{
	     vec[deg_y].resize(deg_x);
	     vec[deg_y].push_back(g) ;
	   }
	   
	 }
	 f_node=rs_export_list_mon_pp_bz_nextnode(f_node);
	 
       }
     return vec;
   }

 
 
//this function store the resulting RURs from RS resolution in a suitable data structure.the latter is used afterward to construct the set of the rurs

std::vector< std::vector< std::vector<Gmpz> > > Rurs_sys_list()
{
  
  std::vector< std::vector< std::vector<Gmpz> > > vect;
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

//this function compute the 2d boxes corresponding to the solutions of a rur
 template< class OutputIterator>
 void affiche_sols_eqs(OutputIterator sol)
     /*     cette fonction fonctionne donc pour un polynome en une variable
            autant que pour une système d'équations*/
 {
   CGALRS_PTR(ident_sols_eqs);
   CGALRS_PTR(ident_node);
   CGALRS_PTR(ident_vect);
   CGALRS_PTR(ident_elt_x);
   CGALRS_PTR(ident_elt_y);
   
   int i,nb_elts;
   /* les solutions des systèmes d'équations (ou des polynômes en une
      variables calculées par RS sont toujours dans une variable par
      defaut : */
   ident_sols_eqs=rs_get_default_sols_eqs();
  /* les solutions sont une liste de vecteurs de longueur le nombre de
     variables, les entrées de ces vecteurs sont des intervalles
     de type mpfi*/
  /* on attrappe le nombre d'éléments de la liste de RS */
  nb_elts=rs_export_list_vect_ibfr_nb(ident_sols_eqs);
  /* on attrappe le premier element de la liste de RS */
  ident_node=rs_export_list_vect_ibfr_firstnode(ident_sols_eqs);
  for (i=1;i<nb_elts+1;i++) {
    /* on dépacte le i-eme élément de la liste */
    ident_vect=rs_export_list_vect_ibfr_monnode(ident_node);
    ident_elt_x=rs_export_elt_vect_ibfr(ident_vect,0);
    ident_elt_y=rs_export_elt_vect_ibfr(ident_vect,1);    
    
    *sol++ = std::make_pair(CGAL::Gmpfi((mpfi_ptr)rs_export_ibfr_mpfi(ident_elt_x)),
			    CGAL::Gmpfi((mpfi_ptr)rs_export_ibfr_mpfi(ident_elt_y)));
    // affiche_vect_ibfr(ident_vect);
    /* on passe a l'élément suivant */
    ident_node=rs_export_list_vect_ibfr_nextnode(ident_node);
  }
 }
 

  } // namespace RS3
 
} // namespace CGAL

#endif  // CGAL_RS_RS3_CALLS_H
