#pragma once

#include "Vector.h"
#include <functional>

namespace kana
{
	template<typename T>
	class Matrix
	{
	public:
		Matrix();
		Matrix(const size_t num_rows, const size_t num_cols, const T val);
		Matrix(const size_t num_rows, const size_t num_cols, const std::function<T()> f);
		virtual ~Matrix() = default;

		Matrix<T>(const Matrix<T>& other) = default;
		Matrix<T>& operator=(const Matrix<T>& other) = default;

		Matrix<T>(Matrix<T>&& other) noexcept;
		Matrix<T>& operator=(Matrix<T>&& other) noexcept;

		Matrix<T>& operator-();

		Matrix<T> operator+(const Matrix<T>& other) const;
		Matrix<T>& operator+=(const Matrix<T>& other);

		Matrix<T> operator-(const Matrix<T>& other) const;
		Matrix<T>& operator-=(const Matrix<T>& other);

		Matrix<T> operator*(const T val) const;
		Matrix<T>& operator*=(const T val);

		Matrix<T> operator/(const T val) const;
		Matrix<T>& operator/=(const T val);

		std::pair<size_t, size_t> shape() const;

		T& at(const size_t row, const size_t col);
		T at(const size_t row, const size_t col) const;
		const Vector<T>& get_row(const size_t ind) const;
		Vector<T> get_col(const size_t ind) const;
		Matrix<T> copy_values() const;

		void print() const;
	private:
		Vector<Vector<T>> m_data;
	};
	
	template<typename T>
	Matrix<T> identity_matrix(size_t size);
}