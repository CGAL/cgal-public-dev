#pragma once
#include "defs.h"
#include <opencv2/core/core.hpp>


using cv::Vec3b;

namespace ColorMap 
{
    /*inline void rainbow(double x, CGAL_Color & color)
	{
		// This function returns a color that varies from blue (0) to red (1)
		// using the usual color spectrum
        double y = 2 * jclamp(0, x, 1) - 1;
		double r = 0, g = 0, b = 0;
		if (y < -0.75) {
			b = 2.5 + 2.0 * y;
		} else if (y < -0.25) {
			b = 1.0;
			g = 1.5 + 2.0 * y;
		} else if (y < 0.25) {
			b = 0.5 - 2.0 * y;
			g = 1.0;
			r = 2.0 * x + 0.5;
		} else if (y < 0.75) {
			g = 1.5 - 2.0 * y;
			r = 1.0;
		} else {
			r = 2.5 - 2.0 * y;
		}
		color = CGAL_Color(uchar(255 * r), uchar(255 * g), uchar(255 * b));
    }*/


	inline Vec3b rainbow(double x)
	{
        double y = 2 * jclamp(0, x, 1) - 1;
		double r, g, b;
		r = g = b = 0;
		if (y < -0.75) {
			b = 2.5 + 2.0 * y;
		} else if (y < -0.25) {
			b = 1.0;
			g = 1.5 + 2.0 * y;
		} else if (y < 0.25) {
			b = 0.5 - 2.0 * y;
			g = 1.0;
			r = 2.0 * y + 0.5;
		} else if (y < 0.75) {
			g = 1.5 - 2.0 * y;
			r = 1.0;
		} else {
			r = 2.5 - 2.0 * y;
		}
		return Vec3b(uchar(255 * b), uchar(255 * g), uchar(255 * r));
	}


    /*inline void blue_white_red(double x, CGAL_Color & color)
	{
		uchar r, g, b;
        double y = jclamp(0, x, 1);
		if (y < 0.5) {
			b = 255;
			g = uchar(510 * y);
			r = uchar(510 * y);
		} else {
			b = uchar(510 * (1 - y));
			g = uchar(510 * (1 - y));
			r = 255;
		}
		color = CGAL_Color(r, g, b);
    }*/


	inline Vec3b blue_white_red(double x)
	{
		uchar r, g, b;
        double y = jclamp(0, x, 1);
		if (y < 0.5) {
			b = 255;
			g = uchar(510 * y);
			r = uchar(510 * y);
		} else {
			b = uchar(510 * (1 - y));
			g = uchar(510 * (1 - y));
			r = 255;
		}
		return Vec3b(r, g, b);
	}
};
