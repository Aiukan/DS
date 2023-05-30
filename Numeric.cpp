#pragma once

#include "Numeric.h"

namespace kana
{
	double differential(std::function<double(double)> f, double point)
	{
		return (f(point + EPSILON) - f(point)) / EPSILON;
	}

	double differential(std::function<std::complex<double>(std::complex<double>)> f, double point)
	{
		return f(std::complex<double>(point, EPSILON)).imag() / EPSILON;
	}
}