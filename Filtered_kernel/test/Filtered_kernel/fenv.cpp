
#include <iostream>

//  https://docs.microsoft.com/en-us/cpp/build/reference/floating-point-optimization
// explains that the below code without fenv_access (on) ignores rhe control_fp
// and might result in cUpper = cLower   (what I cannot reproduce)

// fenv_access (on) requires fp:precise (or fp:strict) 
// that is without the float_control this does not compile with fp:fast
// and we get
// error C3198: invalid use of floating-point pragmas: fenv_access pragma operates only in precise mode

// https://msdn.microsoft.com/en-us/library/45ec64h6.aspx
// explains how to change from fast to precise

namespace Interval {

#pragma float_control(precise, on, push) 

  // With fp:strict fenv is already on
  // With fp:precise we have to do it
#pragma fenv_access (on)  // not in the global scope but in this namespace

  void fct()
{
  std::cout.precision(17);
  double a, b, cLower, cUpper;

  std::cin >> a >> b;

  {
  _controlfp( _RC_DOWN, _MCW_RC );    // round to -infinity
  cLower = a*b;
  _controlfp( _RC_UP, _MCW_RC );       // round to +infinity
  cUpper = a*b;
  _controlfp( _RC_NEAR, _MCW_RC );    // restore rounding mode
  }
  std::cout << "cLower = " << cLower << std::endl;
  std::cout << "cUpper = " << cUpper << std::endl;
}

#pragma float_control(pop)
  
}  // namespace Interval

  
int main()
{
  Interval::fct();
  
 return 0;
}
