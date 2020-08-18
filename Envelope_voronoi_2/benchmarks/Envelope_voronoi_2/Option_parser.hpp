#ifndef OPTION_PARSER_HPP
#define OPTION_PARSER_HPP

#include <string>
#include <vector>
#include <list>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>

#include <CGAL/Benchmark/Option_parser.hpp>

#define DEF_WIN_WIDTH   512
#define DEF_WIN_HEIGHT  512

namespace cb = CGAL::benchmark;
namespace po = boost::program_options;
namespace fi = boost::filesystem;

class Option_parser : public cb::Option_parser 
{
  static const char *ENVELOPE_VERSION;
public:
  /*! Type code */
  enum Type_code {
    TYPE_ENVELOPE_VORONOI = 0,
    TYPE_TRIANGULATION_VORONOI,
    TYPE_ENVELOPE_APOLLONIUS,
    TYPE_CGAL_APOLLONIUS,
    TYPE_EXACUS_APOLLONIUS,
    TYPE_EXACUS_QDX_APOLLONIUS,
    TYPE_SPHERE,
    TYPE_DISPLAY,
    TYPE_E,
    TYPE_T,
    TYPE_A,
    TYPE_C,
    TYPE_X,
    TYPE_Q,
    TYPE_S,
    TYPE_D
  };

  template <class T> class Id {
  public:
    Id() : m_id(0) {}
    Id(unsigned int id) : m_id(id) {}
    unsigned int m_id;
  };

  typedef Id<Type_code>                 Type_id;
  typedef std::vector<Type_id>          Vector_type_id;
  typedef Vector_type_id::iterator      Vector_type_id_iter;

public:
  /*! \brief obtains number of type options */
  static unsigned int get_number_opts(Type_id &);

  /*! Compare the i-th type option to a given option */
  static bool compare_opt(unsigned int i, const char * opt, Type_id &)
  { return strcmp(s_type_opts[i], opt) == 0; }
  
  Option_parser();
  
  enum Generic_option_id { AUTHOR, HELP, LICENSE, VERSION };
  struct Generic_option_exception {
    Generic_option_exception(Generic_option_id option) : m_option(option) {}
    Generic_option_id m_option;
  };

  enum Error_id { FILE_NOT_FOUND, FILE_CANNOT_OPEN };
  struct Error_exception {
    Error_exception(Error_id err) : m_error(err) {}
    Error_id m_error;
  };

  struct Input_file_missing_error : public po::error {
    Input_file_missing_error(std::string & str) : error(str) {}
  };
  
  /*! Parse the options */
  void operator()(int argc, char * argv[]);
  
  /*! Obtain the verbosity level */
  unsigned int get_verbose_level() const { return m_verbose_level; }

  /*! Obtain the number of input files */
  unsigned int get_number_files() const { return m_number_files; }
  
  /*! \brief obtains the base file-name */
  const std::string & get_file_name(unsigned int i) const;

  /*! \brief obtains the full file-name */
  const std::string & get_full_name(unsigned int i) const;

  bool get_postscript() const { return m_postscript; }
  unsigned int get_type_mask() const { return m_type_mask; }

  const char * get_type_name(Type_code id) const { return s_type_opts[id]; }

  /*! Obtain the window width */
  unsigned int get_width() const { return m_win_width; }

  /*! Obtain the window height */
  unsigned int get_height() const { return m_win_height; }
  
  template <class MyId>
  static void my_validate(boost::any & v,
                          const std::vector<std::string> & values);
  
protected:
  /*! The variable map */
  po::variables_map m_variable_map;

  /*! The generic options */
  po::options_description m_generic_opts;

  /*! Visible options */
  po::options_description m_visible_opts;
  
  /*! Command line options */
  po::options_description m_cmd_line_opts;

  /*! Config file options */
  po::options_description m_config_file_opts;

  /*! Environment variable options */
  po::options_description m_environment_opts;

  /*! The configuration option description */
  po::options_description m_config_opts;

  /*! The hidden option description */
  po::options_description m_hidden_opts;    

  /*! Positional option description */
  po::positional_options_description m_positional_opts;

private:
  typedef std::list<fi::path>           Path_list;
  typedef Path_list::iterator           Path_iter;

  typedef std::vector<std::string>      Input_path;
  typedef Input_path::const_iterator    Input_path_const_iterator;
  
  Input_path_const_iterator dirs_begin()
  { return m_variable_map["input-path"].as<Input_path>().begin(); }

  Input_path_const_iterator dirs_end()
  { return m_variable_map["input-path"].as<Input_path>().end(); }

  template <class UnaryFunction>
  UnaryFunction for_each_dir(UnaryFunction func)
  {
    if (!m_variable_map.count("input-path")) return func;
    return std::for_each(dirs_begin(), dirs_end(), func);
  }
  
  /*! A functor that adds a directory to the directory-search structure */
  struct Add_dir {
    Path_list & m_dirs;
    Add_dir(Path_list & dirs) : m_dirs(dirs) {}
    void operator()(const std::string & dir) { m_dirs.push_back(dir); }
  };

  /*! Type options */
  static char * s_type_opts[];

  /*! Verbosity level */
  unsigned int m_verbose_level;  

  /*! The window width */
  unsigned int m_win_width;

  /*! The window height */
  unsigned int m_win_height;

  unsigned int m_type_mask;
  
  bool m_postscript;

  unsigned int m_number_files;

  /*! A collection of directories to search files in */
  Path_list m_dirs;

  std::vector<std::string> m_full_names;
};

#endif