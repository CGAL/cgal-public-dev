/* Only for internal debug purposes may be removed later. */
#define CGAL_PRINT_COLOR 1

#if CGAL_PRINT_COLOR
#define CGAL_BOLD       1
#define CGAL_DARK       2
#define CGAL_UNDERLINE  4
#define CGAL_BLINK      5
#define CGAL_NEGATIVE   7
#define CGAL_BLACK      30
#define CGAL_RED        31
#define CGAL_GREEN      32
#define CGAL_YELLOW     33
#define CGAL_BLUE       34
#define CGAL_MAGNETA    35
#define CGAL_CYAN       36
#define CGAL_WHITE      37

inline std::string change_color_rec(std::ostringstream& o, int color)
{
  o << "\033[1;" << color << "m" << "" << "\033[0m";
  return o.str();
}

template <typename T, typename... Args>
std::string change_color_rec(std::ostringstream& o, int color, const T& value, 
                             Args&... args)
{
  o << "\033[1;" << color << "m" << value << "\033[0m";
  return change_color_rec(o, color, args...);
}

template<typename T, typename... Args>
std::string change_color(int color, const T& value,const Args&... args)
{
  std::ostringstream o;
  o << "\033[1;" << color << "m" << value << "\033[0m";
  return change_color_rec(o, color, args...);
}


#endif
