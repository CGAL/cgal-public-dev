#ifndef CGAL_SHAPE_REGULARIZATION_SAVER_SEGMENTS_2
#define CGAL_SHAPE_REGULARIZATION_SAVER_SEGMENTS_2

// #include <CGAL/license/Shape_regularization.h>

#include <fstream>
#include <vector>
#include <string>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits>
  struct Saver_segments_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Segment = typename GeomTraits::Segment_2;

    inline std::string data() const {
      return out.str();
    }

    void save_segments(const std::vector<Segment> & segments, const std::string &file_name) {
        clear();

        size_t size = 0;
        for (const auto & segment : segments) {
          out << "v " << segment.source() << " " << 0 << std::endl;
          out << "v " << segment.target() << " " << 0 << std::endl;
          out << "v " << segment.target() << " " << 0 << std::endl;
          ++size;
        }

        for (size_t i = 0; i < size * 3; i += 3)
          out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << std::endl;
        
        save(file_name, ".obj");
    }

  private:
    std::stringstream out;

    void clear() {
      out.str(std::string());
    }

    void save(const std::string &file_name, const std::string &extension = ".log") const {
      const std::string final_path = file_name + extension;
      std::ofstream file(final_path.c_str(), std::ios_base::out);

      if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << file_name << std::endl << std::endl;

      file << data() << std::endl;
      file.close();
    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_SAVER_SEGMENTS_2
