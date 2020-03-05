#pragma once
#include "defs_cgal.h"

namespace Skippy {

	inline FT approximate_sqrt(const FT & x)
	{
		// In this function we aim at determining an approximation of S = sqrt(x),
		// that's why we apply the following algorithm :
		// 1. Identification of y0 such that y0^2 ~ S.
		// 2. While yn hasn't converged, y_{n+1} = 0.5 (yn + S / yn).

		// Part 1.
		// As input we have a rational, x = n / d.
		// We get its integral part, s and compute an approximation y0 of sqrt(s).

		CGAL::Gmpz n = x.exact().numerator();
		CGAL::Gmpz d = x.exact().denominator();
		CGAL::Gmpz s = n / d;

		// Gets the approximate number of digits used to represent s.
		// If s ~ 10^m, sqrt(s) ~ 10^(m/2).

		size_t digits = s.approximate_decimal_length();
		std::string str(digits / 2 + 1, '0');
		str[0] = '1';

		FT y_0(str);

		// Part 2.
		// Repeat until convergence

		FT eps_required = FT(1) / FT(10000000000);
		FT eps_curr = 1;
		FT y_curr = y_0, y_next;

		while (eps_curr > eps_required) {
			y_next = (y_curr + x / y_curr) / 2;
			eps_curr = CGAL::abs(y_next - y_curr);
			y_curr = y_next;
		}

		return y_next;
	}
}