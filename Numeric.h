#pragma once

#include <complex>
#include <functional>

namespace kana 
{
	const double PI = 3.14159265359;

	const double EPSILON = 1e-9;

	double differential(std::function<double(double)> f, double point);

	double differential(std::function<std::complex<double>(std::complex<double>)> f, double point);
}