#pragma once

#include <functional>

#include "Vector.h"

namespace kana
{
	std::pair<double, double> bracket_minimum(const std::function<double(double)> f, const double initial = 0, const double s = 1e-2, const double k = 2);

	std::pair<double, double> golden_section_search(const std::function<double(double)> f, const std::pair<double, double>& interval, const size_t n);

	std::pair<double, double> quadratic_fit_search(const std::function<double(double)> f, const double a, const double b, const double c, const size_t n);

	std::pair<double, double> _get_sp_intersection(const double a, const double y_a, const double b, const double y_b, const double l);
	
	std::pair<double, std::vector<std::pair<double, double>>> shubert_piyavskii(const std::function<double(double)> f, const std::pair<double, double>& interval, const double l, const double eps, const double delta = 0.01);

	std::pair<double, double> bisection(const std::function<double(double)> f, const std::pair<double, double>& interval, const double eps = 0.01);
}