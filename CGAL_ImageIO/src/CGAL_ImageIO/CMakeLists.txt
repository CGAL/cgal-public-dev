include(CGAL_library)

CGAL_define_library(
  NAME CGAL_ImageIO
  BUILD_DEFAULT ON
  REQUIRED_DEPENDENCIES CGAL OpenGL
  OPTIONAL_DEPENDENCIES VTK ZLIB
  SOURCES analyze.cpp analyze.h bmp.cpp bmpendian.cpp bmpendian.h bmp.h 
          bmpread.cpp bmpread.h bmptypes.h CMakeLists.txt convert.cpp convert.h 
          fgetns.cpp fgetns.h gif.cpp gif.h gis.cpp gis.h Image_3.cpp ImageIO.cpp 
          inr.cpp inr.h iris.cpp iris.h LICENSE mincio.cpp mincio.h pnm.cpp pnm.h 
          recbuffer.cpp recbuffer.h recline.cpp recline.h reech4x4.cpp reech4x4.h typedefs.h)
