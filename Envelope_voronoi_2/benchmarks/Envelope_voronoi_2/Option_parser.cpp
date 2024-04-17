#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>

#include "Option_parser.hpp"

const char * Option_parser::ENVELOPE_VERSION = "1.0";
char * Option_parser::s_type_opts[] = {
  "envelope_voronoi", "triangulation_voronoi", "envelope_apollonius", 
  "cgal_apollonius", "EXACUS_apollonius", "EXACUS_Qdx_apollonius", "Sphere", 
  "display", "e", "t", "a", "c", "X", "Q", "S", "d"
};

template <class MyId>
void Option_parser::my_validate(boost::any & v,
                                const std::vector<std::string> & values)
{
  typedef std::vector<MyId>     Vector_id;
  MyId tag;
  for (unsigned int i = 0; i < get_number_opts(tag); ++i) {
    if (compare_opt(i, values[0].c_str(), tag)) {
      if (v.empty()) {
        Vector_id vec;
        vec.push_back(MyId(i));
        v = boost::any(vec);
      } else {
        Vector_id vec = boost::any_cast<Vector_id>(v);
        vec.push_back(MyId(i));
        v = boost::any(vec);
      }
      return;
    }
  }
  throw po::validation_error("invalid value");
}

/* Overload the 'validate' function for the user-defined class */
void validate(boost::any & v, const std::vector<std::string> & values,
              Option_parser::Vector_type_id * target_type, int)
{
  Option_parser::my_validate<Option_parser::Type_id>(v, values);
}

/*! Constructor */
Option_parser::Option_parser() :
  m_generic_opts("Generic options"),
  m_config_opts("Configuration options"),
  m_hidden_opts("Hidden options"),
  m_verbose_level(0),
  m_win_width(DEF_WIN_WIDTH),
  m_win_height(DEF_WIN_HEIGHT),
  m_type_mask(0xf),
  m_postscript(false),
  m_number_files(0)
{
  m_dirs.push_front(".");

  const char * root_c = getenv("ROOT");
  if (root_c) 
  {
    fi::path root(root_c, fi::native);
    fi::path dir;
    dir = root / "data";
    m_dirs.push_back(dir);
  }

  // Generic options:
  m_generic_opts.add_options()
    ("author,a", "print author name(s)")
    ("help,h", "print help message")
    ("license,l", "print licence information")
    ("version,v", "print version string")
    ;
    
  // Options allowed on the command line, config file, or env. variables
  m_config_opts.add_options()
    ("input-path,P", po::value<Input_path>()->composing(), "input path")
    ("verbose,V", po::value<unsigned int>(&m_verbose_level)->default_value(0),
     "verbose level")
    ("type", po::value<std::vector<Type_id> >()->composing(),
     "Type\n"
     "  e[nvelope_voronoi]            (0x1)\n"
     "  t[riangulation_voronoi]       (0x2)\n"
     "  [envelope_]a[pollonius]       (0x4)\n"
     "  c[gal_apollonius]             (0x8)\n"
     "  [E]X[ACUS_apollonius]         (0x10)\n"
     "  [EXACUS_]Q[dx_apollonius]     (0x20)\n"
     "  S[phere] - Voronoi on sphere  (0x40)\n"
     "  d[isplay]                     (0x80)\n"
     )
    ("type-mask,T", po::value<unsigned int>(&m_type_mask)->default_value(0x3f),
     "type mask")
    ("width,W",
     po::value<unsigned int>(&m_win_width)->default_value(DEF_WIN_WIDTH ),
     "window width")
    ("height,H",
     po::value<unsigned int>(&m_win_height)->default_value(DEF_WIN_HEIGHT),
     "window height")
    ;
  
  // Options hidden to the user. Allowed only on the command line:
  m_hidden_opts.add_options()
    ("input-file", po::value<Input_path>()->composing(), "input file")
    ;

  m_visible_opts.add(m_generic_opts).add(m_bench_opts).add(m_config_opts);
  m_cmd_line_opts.add(m_generic_opts).add(m_bench_opts).add(m_config_opts).add(m_hidden_opts);
  m_config_file_opts.add(m_bench_opts).add(m_config_opts);
  m_environment_opts.add(m_bench_opts).add(m_config_opts);

  m_positional_opts.add("input-file", -1);  
}

/*! Parse the options */
void Option_parser::operator()(int argc, char * argv[])
{
  po::store(po::command_line_parser(argc, argv).
            options(m_cmd_line_opts).positional(m_positional_opts).run(),
            m_variable_map);
  
  std::ifstream ifs(".bench.cfg");
  po::store(parse_config_file(ifs, m_config_file_opts), m_variable_map);
  po::notify(m_variable_map);

  if (m_variable_map.count("help")) 
  {
    std::cout << m_visible_opts << std::endl;
    throw Generic_option_exception(HELP);
    return;
  }

  if (m_variable_map.count("version")) 
  {
    std::cout << "Voronoi benchmarks v" << ENVELOPE_VERSION << std::endl;
    throw Generic_option_exception(VERSION);
    return;
  }

  if (m_variable_map.count("author")) 
  {
    std::cout << "author Ophir Setter based on code by Efi Fogel" << std::endl;
    throw Generic_option_exception(AUTHOR);
    return;
  }

  if (m_variable_map.count("license"))
  {
    std::cout << "license ???" << std::endl;
    throw Generic_option_exception(LICENSE);
    return;
  }

  if (m_variable_map.count("type"))
  {
    m_type_mask = 0;
    Vector_type_id types = m_variable_map["type"].as<Vector_type_id>();
    for (Vector_type_id_iter it = types.begin(); it != types.end(); ++it) \
    {
      unsigned int size = sizeof(Option_parser::s_type_opts) / 8;
      m_type_mask |= 0x1 << ((*it).m_id % size);
    }
  }
  
  // Add directories specified on the command line via the "input-path" opt.:
  Add_dir tmp = for_each_dir(Add_dir(m_dirs));
  
  if (!m_variable_map.count("input-file")) 
  {
    std::string str("input file missing!");
    throw Input_file_missing_error(str);
    return;
  }

  Input_path files = m_variable_map["input-file"].as<Input_path>();
  m_full_names.resize(files.size());
  Input_path_const_iterator it;
  for (it = files.begin(); it != files.end(); ++it) 
  {
    fi::path file_path(*it, fi::native);
    if (file_path.is_complete())
    {
      if (fi::exists(file_path))
        m_full_names[m_number_files] = file_path.native_file_string();
    }
    else
    {
      for (Path_iter pi = m_dirs.begin(); pi != m_dirs.end(); ++pi)
      {
        fi::path full_file_path = *pi / file_path;
        if (!fi::exists(full_file_path)) continue;
        m_full_names[m_number_files] = full_file_path.native_file_string();
        break;
      }
    }
    
    if (m_full_names[m_number_files].empty())
    {
      std::cerr << "cannot find file " << (*it).c_str() << "!" << std::endl;
      throw Error_exception(FILE_NOT_FOUND);
      return;
    }
    ++m_number_files;
  }
}

/*! Obtain the base file-name */
const std::string & Option_parser::get_file_name(unsigned int i) const
{
  return m_variable_map["input-file"].as<Input_path>()[i];
}

/*! Obtain the full file-name */
const std::string & Option_parser::get_full_name(unsigned int i) const
{ return m_full_names[i]; }

/*! Obtain number of type options */
unsigned int Option_parser::get_number_opts(Type_id &)
{ return sizeof(s_type_opts) / sizeof(char *); }