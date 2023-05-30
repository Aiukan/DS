#pragma once

#include <functional>
#include "Vector.h"
#include "Matrix.h"

namespace kana
{
	template <typename T>
	T min(const Vector<T>& vec);

	template <typename T>
	T max(const Vector<T>& vec);

	template <typename T>
	T arithmetic_mean(const Vector<T>& vec);

	template <typename T>
	T dot(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	T magnitude(const Vector<T>& vec);

	template <typename T>
	T distance(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	T median(const Vector<T>& vec);

	template <typename T>
	T quantile(const Vector<T>& vec, const double p);

	template <typename T>
	T mode(const Vector<T>& vec);

	template <typename T>
	T range(const Vector<T>& vec);

	template <typename T>
	Vector<T> demeaning_mean(const Vector<T>& vec);

	template <typename T>
	T variance(const Vector<T>& vec);

	template <typename T>
	T standard_deviation(const Vector<T>& vec);

	template <typename T>
	T interquartile_range(const Vector<T>& vec);

	template <typename T>
	T covariance(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	T correlation(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	Vector<T> bucketize(const Vector<T>& vec, const size_t bucket_num);

	template <typename T>
	T directional_variance(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	Vector<T> directional_variance_gradient(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	Vector<T> direction(const Vector<T>& vec);

	template <typename T>
	Vector<T> project(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	Vector<T> remove_projection(const Vector<T>& vec1, const Vector<T>& vec2);

	template <typename T>
	Vector<T> fold(const Matrix<T>& mat);

	template <typename T>
	Vector<T> fold_mean(const Matrix<T>& mat);

	template <typename T>
	Matrix<T> correlation_matrix(const Matrix<T>& mat);

	template <typename T>
	bool is_identity(const Matrix<T>& mat);

	template <typename T>
	Matrix<T> de_meaned(const Matrix<T>& mat);

	template <typename T>
	Matrix<T> rescaled(const Matrix<T>& mat);

	template <typename T>
	T directional_variance(const Matrix<T>& mat, const Vector<T>& vec);

	template <typename T>
	Vector<T> directional_variance_gradient(const Matrix<T>& mat, const Vector<T>& vec);

	template <typename T>
	Vector<T> first_principal_component(const Matrix<T>& mat);

	template <typename T>
	Matrix<T> remove_projection(const Matrix<T>& mat, const Vector<T>& vec);

	double uniform_pdf(const double x);

	double uniform_cdf(const double x);

	double normal_pdf(const double x, const double mu = 0, const double sigma = 1);

	double normal_cdf(const double x, const double mu = 0, const double sigma = 1);

	double inverse_normal_cdf(const double p, const double mu = 0,
		const double sigma = 1, const double tolerance = 0.00001);

	int bernoulli_trial(const double p);

	int binomial(const int n, const double p);

	double normal_lower_bound(const double probability, const double mu = 0, const double sigma = 1);

	double normal_upper_bound(const double probability, const double mu = 0, const double sigma = 1);

	std::pair<double, double> normal_bounds(const double probability, const double mu = 0, const double sigma = 1);

	double beta_pdf(const double x, const double alpha, const double beta);

	template <typename T>
	Vector<T> gradient_descent(const std::function<T(const Vector<T>&)> f_target, const std::function<Vector<T>(Vector<T>)> f_gradient,
		const Vector<T>& starting_point, bool maximize = false, const double tolerance = 0.000001);

	template <typename T>
	Vector<T> gradient_descent_stochastic(const std::function<T(Vector<T>, Vector<T>, Vector<T>)> f_target,
		const std::function<Vector<T>(Vector<T>, Vector<T>, Vector<T>)> f_gradient,
		Matrix<T>& x, Matrix<T>& y, const Vector<T>& starting_point, bool maximize = false, const double alpha_0 = 0.01);
}