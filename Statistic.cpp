#pragma once

#include "Statistic.h"
#include "Numeric.h"
#include "Random.h"
#include <cmath>

namespace kana 
{
	template <typename T>
	T min(const Vector<T>& vec)
	{
		return *std::min_element(cbegin(vec), cend(vec));
	}

	template <typename T>
	T max(const Vector<T>& vec)
	{
		return *std::max_element(cbegin(vec), cend(vec));
	}

	template <typename T>
	T arithmetic_mean(const Vector<T>& vec)
	{
		T sum = 0;
		for (size_t i = 0; i < vec.size(); ++i)
			sum += vec[i];
		return vec.size();
	}

	template <typename T>
	T dot(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		T sum = 0;
		for (size_t i = 0; i < vec1.size(); ++i)
			sum += vec1[i] * vec2[i];
		return sum;
	}

	template <typename T>
	T magnitude(const Vector<T>& vec)
	{
		return std::sqrt(dot(vec, vec));
	}

	template <typename T>
	T distance(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		Vector<T> diff = vec1 - vec2;
		return diff.magnitude();
	}

	template <typename T>
	T median(const Vector<T>& vec)
	{
		Vector<T> new_vec = vec.sorted();
		if (new_vec.size() % 2 == 1)
			return new_vec[new_vec.size() / 2];
		else
			return (new_vec[new_vec.size() / 2 - 1] + new_vec[new_vec.size() / 2]) / 2;
	}

	template <typename T>
	T quantile(const Vector<T>& vec, const double p)
	{
		Vector<T> new_vec = vec.sorted();
		if (p == 1.0)
			return new_vec[new_vec.size() - 1];
		size_t ind = static_cast<size_t>(p * new_vec.size());
		return new_vec[ind];
	}

	template <typename T>
	T mode(const Vector<T>& vec)
	{
		Vector<T> new_vec = vec.sorted();
		T val = new_vec[0];
		size_t k = 1, maxk = 1;
		for (size_t i = 1; i < new_vec.size(); ++i) {
			if (new_vec[i] == new_vec[i - 1])
				++k;
			else {
				if (k > maxk) {
					maxk = k;
					val = new_vec[i - 1];
				}
				k = 1;
			}
		}
		if (k > maxk) {
			val = new_vec[new_vec.size() - 1];
		}
		return val;
	}

	template <typename T>
	T range(const Vector<T>& vec)
	{
		T min_val = vec.min();
		T max_val = vec.max();
		return max_val - min_val;
	}

	template <typename T>
	Vector<T> demeaning_mean(const Vector<T>& vec)
	{
		T mean = vec.arithmetic_mean();
		Vector<T> mean_vec(vec.size(), mean);
		Vector<T> new_vec = vec - mean_vec;
		return new_vec;
	}

	template <typename T>
	T variance(const Vector<T>& vec)
	{
		Vector<T> de_vec = vec.demeaning_mean();
		return de_vec.dot(de_vec) / (de_vec.size() - 1);
	}

	template <typename T>
	T standard_deviation(const Vector<T>& vec)
	{
		return std::sqrt(vec.variance());
	}

	template <typename T>
	T interquartile_range(const Vector<T>& vec)
	{
		return vec.quantile(0.75) - vec.quantile(0.25);
	}

	template <typename T>
	T covariance(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		Vector<T> de_vec1 = vec1.demeaning_mean();
		Vector<T> de_vec2 = vec2.demeaning_mean();
		return de_vec1.dot(de_vec2) / (de_vec1.size() - 1);
	}

	template <typename T>
	T correlation(const Vector<T>& vec1, const Vector<T> vec2)
	{
		T std_dev1 = vec1.standard_deviation();
		T std_dev2 = vec2.standard_deviation();
		if (std_dev1 > 0 && std_dev2 > 0)
			return covariance(vec1, vec2) / std_dev1 / std_dev2;
		return 0;
	}

	template <typename T>
	Vector<T> bucketize(const Vector<T>& vec, const size_t bucket_num)
	{
		double min_val = vec.min();
		double max_val = vec.max();
		if (min_val == max_val)
			return Vector(1, vec.size());
		double interval = (max_val - min_val) / bucket_num;
		size_t ind;
		Vector<T> bucketed = Vector<T>(bucket_num, 0);
		for (size_t i = 0; i < vec.size(); ++i) {
			ind = static_cast<size_t>((vec[i] - min_val) / interval);
			if (ind == bucket_num)
				ind--;
			bucketed[ind]++;
		}
		return bucketed;

	}

	template <typename T>
	Vector<T> direction(const Vector<T>& vec)
	{
		Vector<T> res = vec.copy_values();
		T len = magnitude(res);
		return res / len;
	}

	template <typename T>
	T directional_variance(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		T direc_len = dot(vec1 * direction(vec2));
		return direc_len * direc_len;
	}

	template <typename T>
	Vector<T> directional_variance_gradient(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		T direc_len = dot(vec1 * direction(vec2));
		Vector<T> gradient = new Vector<T>();
		for (size_t i = 0; i < vec1.size(); ++i)
			gradient.push_back(vec1[i] * (2 * direc_len));
		return gradient;
	}

	template <typename T>
	Vector<T> project(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		T len = dot(vec1, vec2);
		Vector<T> res = vec2.copy_values();
		return res * len;
	}

	template <typename T>
	Vector<T> remove_projection(const Vector<T>& vec1, const Vector<T>& vec2)
	{
		Vector<T> res = vec1.copy_values();
		return res - res.project(vec2);
	}

	template<typename T>
	Vector<T> fold(const Matrix<T>& mat)
	{
		Vector<T> new_vec(mat.get_row(0));
		for (size_t i = 1; i < mat.shape().first; ++i)
			for (size_t j = 0; j < mat.shape().second; ++j)
				new_vec[j] += mat.at(i, j);
		return new_vec;
	}

	template<typename T>
	Vector<T> fold_mean(const Matrix<T>& mat)
	{
		Vector<T> new_vec(mat.get_row(0));
		for (size_t i = 1; i < mat.shape().first; ++i)
			for (size_t j = 0; j < mat.shape().second; ++j)
				new_vec[j] += mat.at(i, j);
		new_vec /= mat.shape().first;
		return new_vec;
	}

	template<typename T>
	bool is_identity(const Matrix<T>& mat)
	{
		auto shape = mat.shape();
		if (shape.first != shape.second)
			return false;
		for (int i = 0; i < shape.first; ++i)
			for (int j = 0; j < shape.first; ++j)
				if (i == j) {
					if (mat.at(i, j) != 1.0)
						return false;
				}
				else if (mat.at(i, j) != 0.0)
					return false;
		return true;
	}

	template<typename T>
	Matrix<T> correlation_matrix(const Matrix<T>& mat)
	{
		Matrix<T> result = Matrix<T>(mat.shape().second, mat.shape().second, 1);
		Vector<T> first;
		Vector<T> second;
		T corr;
		for (size_t i = 0; i < mat.shape().second - 1; ++i)
			for (size_t j = i + 1; j < mat.shape().second; ++j) {
				Vector<T> first = mat.get_col(i);
				Vector<T> second = mat.get_col(j);
				corr = correlation(first, second);
				result.at(i, j) = corr;
				result.at(j, i) = corr;
			}
		return result;
	}

	template<typename T>
	Matrix<T> de_meaned(const Matrix<T>& mat)
	{
		Vector<T> col;
		Matrix<T> result = mat.copy_values();
		T mean;
		for (size_t j = 0; j < result.shape().second; ++j) {
			col = result.get_col(j);
			mean = col.arithmetic_mean();
			for (size_t i = 0; i < result.shape().first; ++i)
				result.at(i, j) -= mean;
		}
		return result;
	}

	template<typename T>
	Matrix<T> rescaled(const Matrix<T>& mat)
	{
		Vector<T> col;
		Matrix<T> result = mat.copy_values();
		T mean;
		T std_deviation;
		size_t n = result.shape().first;
		size_t m = result.shape().second;
		for (size_t j = 0; j < m; ++j) {
			col = result.get_col(j);
			std_deviation = col.standard_deviation();
			if (std_deviation > 0) {
				mean = col.arithmetic_mean();
				for (size_t i = 0; i < n; ++i) {
					result.at(i, j) -= mean;
					result.at(i, j) / std_deviation;
				}
			}
		}
		return result;
	}

	template<typename T>
	T directional_variance(const Matrix<T>& mat, const Vector<T>& vec)
	{
		T len = 0;
		Vector<T> dir = direction<T>(vec);
		for (size_t i = 0; i < mat.shape().first; ++i)
			len += dot(mat.get_row(i), dir);
		return len * len;
	}

	template<typename T>
	Vector<T> directional_variance_gradient(const Matrix<T>& mat, const Vector<T>& vec)
	{
		Vector<T> gradient = Vector<T>(mat.shape().second, 0);
		Vector<T> dir = direction(vec);
		T len;
		for (size_t i = 0; i < mat.shape().first; ++i) {
			len = dot(mat.get_row(i), dir);
			for (size_t j = 0; j < mat.shape().second; ++j)
				gradient[j] += mat.at(i, j) * (2 * len);
		}
		return gradient;
	}

	template<typename T>
	Vector<T> first_principal_component(const Matrix<T>& mat)
	{
		Vector<T> guess = Vector<T>(mat.shape().second, 1);
		auto dir_var = std::bind(kana::directional_variance, mat, std::placeholders::_1);
		auto dir_var_grad = std::bind(directional_variance_gradient, mat, std::placeholders::_1);
		Vector<T> unscaled_maximizer = gradient_descent(dir_var, dir_var_grad, guess);
		return direction(unscaled_maximizer);
	}

	template<typename T>
	Matrix<T> remove_projection(const Matrix<T>& mat, const Vector<T>& vec)
	{
		size_t n = mat.shape().first;
		size_t m = mat.shape().second;
		Matrix<T> res = mat.copy_values();
		Vector<T> row;
		T len;
		for (size_t i = 0; i < n; ++i) {
			len = dot(res.get_row(i), vec);
			for (size_t j = 0; j < m; ++j)
				res.at(i, j) = vec[j] * len;
		}
		return res;
	}


	double uniform_pdf(const double x)
	{
		if (x >= 0 && x < 1)
			return 1;
		return 0;
	}

	double uniform_cdf(const double x)
	{
		if (x < 0)
			return 0;
		if (x < 1)
			return x;
		return 1;
	}

	double normal_pdf(const double x, const double mu, const double sigma)
	{
		double ex = std::exp(-(x - mu) * (x - mu) / 2 / sigma / sigma);
		return ex / std::sqrt(2 * PI) / sigma;
	}

	double normal_cdf(const double x, const double mu, const double sigma)
	{
		return (1 + std::erf((x - mu) / std::sqrt(2) / sigma)) / 2;
	}

	double inverse_normal_cdf(const double p, const double mu,
		const double sigma, const double tolerance)
	{
		if (mu != 0 || sigma != 1)
			return mu + sigma * inverse_normal_cdf(p, 0, 1, tolerance);

		double lower_val = -10.0;
		double upper_val = 10.0;
		double middle_val = 0, middle_prob = 0;
		while (upper_val - lower_val > tolerance)
		{
			middle_val = (upper_val + lower_val) / 2;
			middle_prob = normal_cdf(middle_val);
			if (middle_prob < p)
				lower_val = middle_val;
			else
				upper_val = middle_val;
		}

		return middle_val;
	}

	int bernoulli_trial(const double p)
	{
		if (p >= random_uniform())
			return 1;
		return 0;

	}

	int binomial(const int n, const double p)
	{
		int res = 0;
		for (int i = 0; i < n; ++i)
			res += bernoulli_trial(p);
		return res;
	}

	double normal_lower_bound(const double probability, const double mu, const double sigma)
	{
		return inverse_normal_cdf(1 - probability, mu, sigma);
	}

	double normal_upper_bound(const double probability, const double mu, const double sigma)
	{
		return inverse_normal_cdf(probability, mu, sigma);
	}

	std::pair<double, double> normal_bounds(const double probability, const double mu, const double sigma)
	{
		double lower = inverse_normal_cdf((1 - probability) / 2, mu, sigma);
		double upper = inverse_normal_cdf((probability + 1) / 2, mu, sigma);
		std::pair<double, double> bounds = std::make_pair(lower, upper);
		return bounds;
	}

	double beta_pdf(const double x, const double alpha, const double beta)
	{
		double b = std::tgamma(alpha) * std::tgamma(beta) / (std::tgamma(alpha + beta));
		return std::pow(x, alpha - 1) * std::pow(1 - x, beta - 1) / b;
	}

	template<typename T>
	Vector<T> gradient_descent(const std::function<T(const Vector<T>&)>& f_target, const std::function<Vector<T>(Vector<T>)>& f_gradient,
		const Vector<T>& starting_point, bool maximize, const double tolerance)
	{
		Vector<T> point = starting_point;
		T value = f_target(point);
		Vector<T> gradient;
		Vector<T> next_point;
		T next_value;
		bool changed;
		T change;
		double sign = -1;
		if (maximize)
			sign = 1;

		size_t cycles = 0;
		while (true) {
			gradient = f_gradient(point);
			changed = false;
			if (cycles > 50) {
				cycles = 0;
				for (double i = 1; i < 100000; i *= 2) {
					next_point = point + gradient * (i * sign);
					next_value = f_target(next_point);
					if (next_value * sign < value * sign) {
						if (i != 1)
							changed = true;
						break;
					}
				}
				if (!changed && next_value * sign > value * sign)
					throw std::invalid_argument("Results in infinity");
			}
			double step = 1;
			while (step != 0) {
				next_point = point + gradient * (step * sign);
				next_value = f_target(next_point);
				if (next_value * sign > value * sign)
					break;
				else
					step /= 2;
			}
			change = next_value - value;
			if (change < tolerance && change > -tolerance)
				return point;
			else {
				point = next_point;
				value = next_value;
				cycles++;
			}
		}
		return point;
	}

	template<typename T>
	Vector<T> gradient_descent_stochastic(const std::function<T(Vector<T>, Vector<T>, Vector<T>)> f_target,
		const std::function<Vector<T>(Vector<T>, Vector<T>, Vector<T>)> f_gradient,
		Matrix<T>& x, Matrix<T>& y, const Vector<T>& starting_point, bool maximize, const double alpha_0)
	{
		Vector<size_t> indexes = Vector<size_t>();
		for (int i = 0; i < x.shape().first; ++i)
			indexes.push_back(i);
		Vector<T> point = starting_point;
		T value;
		double sign = -1;
		if (maximize)
			sign = 1;
		Vector<T> best_point = point.copy_values();
		T best_value = std::numeric_limits<T>::infinity() * sign * -1;
		double alpha = alpha_0;
		size_t pointless_cycles = 0;
		while (pointless_cycles < 15) {
			value = 0;
			for (size_t i = 0; i < indexes.size(); ++i)
				value += f_target(x.get_row(i), y.get_row(i), point);
			if (value * sign > best_value * sign) {
				best_point = point;
				best_value = value;
				pointless_cycles = 0;
				alpha = alpha_0;
			}
			else {
				pointless_cycles++;
				alpha *= 0.9;
			}
			shuffle(indexes);
			for (size_t i = 0; i < indexes.size(); ++i) {
				point += f_gradient(x.get_row(indexes[i]), y.get_row(indexes[i]), point) * alpha * sign;
			}
		}
		return point;
	}

}