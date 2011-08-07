#ifndef POLYGON_OFFSET_2_COLORS_H
#define POLYGON_OFFSET_2_COLORS_H

#include <QColor>

struct Polygon_color
{
  QColor polygon;
  QColor kgon;
  QColor kgon_sum;
  QColor kgon_circles;
  QColor kgon_skipped_circles;
  QColor kgon_offset;
  QColor exact_offset;
  QColor approximate_offset;
  QColor skeleton;
  QColor core;
  QColor approximability_positive;
  QColor approximability_negative;
  QColor red;
  QColor green;
  QColor blue;
  QColor yellow;
  QColor cyan;
  QColor magenta;
  QColor cyan_dark;
  QColor magenta_dark;
  QColor grey;
  QColor orange;
  QColor green_dark;
  QColor red_dark;
  
  Polygon_color() :
    red(255,0,0),
    green(0,255,0),
    blue(0,0,255),
    yellow(255,255,0),
    orange(255,127,0),
    cyan(0,255,255),
    magenta(255,0,255),
    cyan_dark(0,127,127),
    magenta_dark(127,0,127),
    grey(180,180,180),
    green_dark(0,127,0),
    red_dark(127,0,0),
    polygon(0,0,255), //blue
    kgon(255,0,0), // red
    kgon_sum(0,255,0), // green
    kgon_circles(1,246,203), // light blue
    kgon_skipped_circles(180,180,180), // grey
    kgon_offset(0,255,127), //(127,255,191),
    exact_offset(255,0,0), // red
    approximate_offset(255,0,255), // magenta
    skeleton(255,0,127), // cold red
    core(255,127,0), // orange
    approximability_positive(64, 64, 255), // bluish
    approximability_negative(191, 191, 0) // yellowish
  {}
/*
  Polygon_color() :
    polygon(0,0,255), // blue
    kgon(0,0,0), // orange
    kgon_sum(255,0,0), // red
    kgon_circles(229,137,133), // light red
    kgon_skipped_circles(180,180,180), // grey
    circle_sum(127,255,191),
    exact_offset(50,153,50), //dark  green
    approximate_offset(85,244,150) // turkis
  {}
*/
};

struct Construction_color
{
  QColor eps;
  QColor rme;
  QColor rad;
  QColor rpe;

  Construction_color() :
  eps(255,181,0),
  rme(251,125,14),
  rad(213,128,52),
  rpe(133,77,28)
  {}
};

static Polygon_color polygon_color;
static Construction_color line_color;
static Construction_color circle_color;

#endif // POLYGON_OFFSET_2_COLORS_H
