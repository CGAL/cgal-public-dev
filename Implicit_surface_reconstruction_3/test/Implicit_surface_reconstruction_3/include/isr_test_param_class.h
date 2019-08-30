// rename file isr_test_parameters.h
// need to be consistent with the delimiter comments + {} style + tabs
// replace BOOST_FOREACH by for( :)
// replace pwnl by point_set or pwn_set or point_set_w_normal
#ifndef ISR_TEST_PARAM_CLASS_H
#define ISR_TEST_PARAM_CLASS_H

//includes
#include <iostream>
#include <list>

//Param struct
struct Param { // rename TestParameter
  bool octree; // later: merge octree and del_ref in indicator_func_refinement as enum type
  bool del_ref;
  bool poisson; // same with equation_to_solve
  bool spectral;
  bool march_tets; // same with contouring_method
  bool make_sm;

  void display(std::ostream &stream) const ;
};

void Param::display(std::ostream &stream) const // in struct directly
{
  bool first = true ;
  stream << "Used features : " ;

  if (octree) 
  {
    if (!first)
      stream << ", " ;
    stream << "octree" ;
    first = false;
  }

  if (del_ref) 
  {
    if (!first)
      stream << ", " ;
    stream << "del_ref" ;
    first = false;
  }

  if (poisson) 
  {
    if (!first)
      stream << ", " ;
    stream << "poisson" ;
    first = false;
  }

  if (spectral) 
  {
    if (!first)
      stream << ", " ;
    stream << "spectral" ;
    first = false;
  }

  if (march_tets) 
  {
    if (!first)
      stream << ", " ;
    stream << "march_tets" ;
    first = false;
  }

  if (make_sm) 
  {
    if (!first)
      stream << ", " ;
    stream << "make_sm" ;
    first = false;
  }
}

std::ostream &operator<<( std::ostream &stream, const Param &p)
{
  p.display(stream);
  return (stream);
}

//Parameters class ----- /!\ parametres avec octree + spectral desactives pour l'instant
class Parameters { // rename TestParameterList

  public :
  //Constructor
  Parameters(bool without_sm = false) // maybe child clase TestParameterListNoCGALMesher
  {
    if(without_sm)
    {
      Param p2 = {false,  true, false, true, true, false};
      Param p4 = {false,  true, true, false, true, false};
      // Param p6 = {true,  false, false, true, true, false};
      Param p8 = {true,  false, true, false, true, false};

      paramList.push_back(p2);
      paramList.push_back(p4);
      // paramList.push_back(p6);
      paramList.push_back(p8);
    }

    else
    {
      Param p1 = {false,  true, false, true, false, true};
      Param p2 = {false,  true, false, true, true, false};
      Param p3 = {false,  true, true, false, false, true};
      Param p4 = {false,  true, true, false, true, false};
      // Param p5 = {true,  false, false, true, false, true};
      // Param p6 = {true,  false, false, true, true, false};
      Param p7 = {true,  false, true, false, false, true};
      Param p8 = {true,  false, true, false, true, false};

      paramList.push_back(p1);
      paramList.push_back(p2);
      paramList.push_back(p3);
      paramList.push_back(p4);
      // paramList.push_back(p5);
      // paramList.push_back(p6);
      paramList.push_back(p7);
      paramList.push_back(p8);
    }
  }

  //public methods
  void add(const Param p)
  {
    paramList.push_back(p);
  }

  std::list<Param>::const_iterator begin() const
  {
    return (paramList.begin());
  }

  std::list<Param>::const_iterator end() const
  {
    return (paramList.end());
  }

  protected :

  //Attributes
  std::list<Param> paramList;
};

//test function
template <typename TestFunctorT>
bool test_all_param(TestFunctorT test_function, PwnList &input_pwnl, bool without_sm = false)
{
  bool success = true;
  bool curr_par_success = true;
  Parameters plist(without_sm);
  for (std::list<Param>::const_iterator param = plist.begin() ; param != plist.end() ; param++) { // for (const TestParameter &param : param_list)
    curr_par_success = true;
    std::cout << "///////////" << " " << *param << " "<< "///////////" << std::endl;
    if (!test_function.run(*param, input_pwnl)) {
      success = false ;
      curr_par_success = false;
    }
    std::cout << "/////////////////////////// " << (curr_par_success ? "PASSED" : "FAILED") << " ///////////////////////////" << std::endl;
    std::cout << std::endl;
  }
  return (success); // return true
}

#endif //ISR_TEST_PARAM_CLASS_H
