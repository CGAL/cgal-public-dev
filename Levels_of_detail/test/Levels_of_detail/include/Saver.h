#ifndef CGAL_LOD_SAVER_H
#define CGAL_LOD_SAVER_H

#if defined(WIN32) || defined(_WIN32)
#define _NL_ "\r\n"
#else
#define _NL_ "\n"
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

namespace CGAL {

  namespace Levels_of_detail {

    class Saver {

    public:
      
      Saver() { 
        out.precision(20); 
      }

      void clear() { 
        out.str(std::string()); 
      }

    private:
      std::stringstream out;

      inline std::string data() const { 
        return out.str(); 
      }

      void save(const std::string file_path) const {
        std::ofstream file(file_path.c_str(), std::ios_base::out);

        if (!file)
          std::cerr << std::endl
                    << "ERROR: Error saving file " << file_path << std::endl
                    << std::endl;

        file << data() << std::endl;
        file.close();
      }
      
    }; // Saver

  } // namespace Levels_of_detail

} // namespace CGAL

#endif // CGAL_LOD_SAVER_H
