#pragma once

#include "Optimization.h"
#include "Numeric.h"
#include <vector>

namespace kana
{
	std::pair<double, double> bracket_minimum(const std::function<double(double)> f, const double initial, const double s, const double k)
	{
		double a = initial, y_a = f(initial);
		double b = initial + s, y_b = f(initial + s);
		double step = s;
		if (y_b > y_a) {
			std::swap(a, b);
			std::swap(y_a, y_b);
			step = -step;
		}
		while (true) {
			double c = b + step, y_c = f(b + step);
			if (y_c > y_b)
				return a < c ? std::make_pair(a, c) : std::make_pair(c, a);
			a = b;
			y_a = y_b;
			b = c;
			y_b = y_b;
			step *= k;
		}
	}

	std::pair<double, double> golden_section_search(const std::function<double(double)> f, const std::pair<double, double>& interval, const size_t n)
	{
		double rho = 1 / PHI;
		double a = interval.first, b = interval.second;
		double d = rho * b + (1 - rho) * a;
		double y_d = f(d);
		for (size_t i = 0; i < n; i++) {
			double c = rho * a + (1 - rho) * b;
			double y_c = f(c);
			if (y_c < y_d) {
				b = d;
				d = c;
				y_d = y_c;
			}
			else {
				a = b;
				b = c;
			}
		}
		return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
	}

	std::pair<double, double> quadratic_fit_search(const std::function<double(double)> f, const double a0, const double b0, const double c0, const size_t n)
	{
		double a = a0, b = b0, c = c0;
		double y_a = f(a), y_b = f(b), y_c = f(c);
		for (size_t i = 0; i < n; i++) {
			double x = 0.5 * (y_a * (b * b - c * c) + y_b * (c * c - a * a) + y_c * (a * a - b * b)) /
				(y_a * (b - c) + y_b * (c - a) + y_c * (a - b));
			double y_x = f(x);
			if (x > b) {
				if (y_x > y_b) {
					c = x;
					y_c = y_x;
				}
				else {
					a = b;
					y_a = y_b;
					b = x;
					y_b = y_x;
				}
			}
			else if (x < b) {
				if (y_x > y_b) {
					a = x;
					y_a = y_x;
				}
				else {
					c = b;
					y_c = y_b;
					b = x;
					y_b = y_x;
				}
			}
		}
		return std::make_pair(a, b);
	}

	std::pair<double, double> _get_sp_intersection(const double a, const double y_a, const double b, const double y_b, const double l)
	{
		double t = ((y_a - y_b) - l * (a - b))/ 2 / l;
		return std::make_pair(a + t, y_a - t * l);
	};
	
	std::pair<double, std::vector<std::pair<double, double>>> shubert_piyavskii(const std::function<double(double)> f, const std::pair<double, double>& interval, const double l, const double eps, const double delta)
	{
		double a = interval.first, y_a = f(a);
		double b = interval.second, y_b = f(b);
		double m = (a + b) / 2, y_m = f(m);
		std::vector<std::pair<double, double>> points;
		points.push_back(std::make_pair(a, y_a));
		points.push_back(_get_sp_intersection(a, y_a, m, y_m, l));
		points.push_back(std::make_pair(m, y_m));
		points.push_back(_get_sp_intersection(m, y_m, b, y_b, l));
		points.push_back(std::make_pair(b, y_b));
		double diff = INFINITY;
		while (diff > eps) {
			double y_min = INFINITY;
			size_t i_min = 0;
			for (size_t i = 0; i < points.size(); i++) 
				if (y_min < points[i].second) {
					i_min = i;
					y_min = points[i].second;
				}
			std::pair<double, double> point_cur = std::make_pair(points[i_min].first, f(points[i_min].first));
			diff = point_cur.second - y_min;
			std::pair<double, double> point_prev = _get_sp_intersection(points[i_min - 1].first, points[i_min - 1].second, point_cur.first, point_cur.second, l);
			std::pair<double, double> point_next = _get_sp_intersection(point_cur.first, point_cur.second, points[i_min + 1].first, points[i_min + 1].second, l);
			points.erase(points.begin() + i_min);
			points.insert(points.begin() + i_min, point_next);
			points.insert(points.begin() + i_min, point_cur);
			points.insert(points.begin() + i_min, point_prev);
		}
		std::vector<std::pair<double, double>> intervals;
		double y_min = INFINITY;
		size_t i_min = 0;
		for (size_t i = 0; i < points.size(); i += 2)
			if (y_min < points[i].second) {
				i_min = i;
				y_min = points[i].second;
			}
		std::pair<double, double> point_min = std::make_pair(points[i_min].first, y_min);
		for (size_t i = 1; i < points.size(); i++) {
			double dy = y_min - points[i].second;
			double x_lo = std::max(a, points[i].first - dy / l);
			double x_hi = std::min(b, points[i].first + dy / l);
			if (intervals.size() != 0 && intervals[intervals.size() - 1].second + delta >= x_lo)
				intervals[intervals.size() - 1].second = x_hi;
			else
				intervals.push_back(std::make_pair(x_lo, x_hi));
		}
		return std::make_pair(y_min, intervals);
	}

	std::pair<double, double> bisection(const std::function<double(double)> f, const std::pair<double, double>& interval, const double eps)
	{
		std::pair<double, double> cur_interval = std::make_pair(interval.first, interval.second);
		if (cur_interval.first >= cur_interval.second)
			std::swap(cur_interval.first, cur_interval.second);
		double y_a = f(cur_interval.first);
		double y_b = f(cur_interval.second);
		if (std::signbit(y_a) == std::signbit(y_b))
			return cur_interval;
		while (cur_interval.second - cur_interval.first > eps) {
			double mid = (cur_interval.first + cur_interval.second) / 2;
			if (std::signbit(y_a) == std::signbit(f(mid))) {
				cur_interval.first = mid;
				y_a = f(cur_interval.first);
			}
			else {
				cur_interval.second = mid;
				y_b = f(cur_interval.second);
			}
		}
		return cur_interval;
	}
}