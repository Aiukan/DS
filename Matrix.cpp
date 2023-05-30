#pragma once

#include "Matrix.h"
#include <iostream>

namespace kana 
{
	template<typename T>
	Matrix<T>::Matrix()
		: m_data()
	{

	}

	template<typename T>
	Matrix<T>::Matrix(const size_t num_rows, const size_t num_cols, const T val)
		: m_data(num_rows, Vector<T>(num_cols, val))
	{
	}

	template<typename T>
	Matrix<T>::Matrix(const size_t num_rows, const size_t num_cols, const std::function<T()> f)
		: m_data(num_rows, Vector<T>(num_cols, f))
	{
	}

	template<typename T>
	Matrix<T>::Matrix(Matrix<T>&& other) noexcept
		: m_data(std::move(other.m_data))
	{
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(Matrix<T>&& other) noexcept
	{
		this->m_data = move(other.m_data);
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator-()
	{
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			m_data[i] = -m_data[i];
		}

		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
	{
		Matrix<T> new_matrix(this->shape().first, this->shape().second, T());
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			new_matrix.m_data[i] = this->m_data[i] + other.m_data[i];
		}
		return new_matrix;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other)
	{
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			this->m_data[i] += other.m_data[i];
		}
		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
	{
		Matrix<T> new_matrix(this->shape().first, this->shape().second, T());
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			new_matrix.m_data[i] = this->m_data[i] - other.m_data[i];
		}
		return new_matrix;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other)
	{
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			this->m_data[i] -= other.m_data[i];
		}
		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator*(const T val) const
	{
		Matrix new_matrix(this->shape().first, this->shape().second, 0);
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			new_matrix.m_data[i] = this->m_data[i] * val;
		}
		return new_matrix;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator*=(const T val)
	{
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			this->m_data[i] *= val;
		}
		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator/(const T val) const
	{
		Matrix new_matrix(this->shape().first, this->shape().second, 0);
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			new_matrix.m_data[i] = this->m_data[i] / val;
		}
		return new_matrix;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator/=(const T val)
	{
		for (size_t i = 0; i < this->shape().first; ++i)
		{
			this->m_data[i] /= val;
		}
		return *this;
	}

	template<typename T>
	std::pair<size_t, size_t> Matrix<T>::shape() const
	{
		return std::make_pair(m_data.size(), m_data[0].size());
	}

	template<typename T>
	T& Matrix<T>::at(const size_t row, const size_t col)
	{
		return m_data[row][col];
	}

	template<typename T>
	T Matrix<T>::at(const size_t row, const size_t col) const
	{
		return m_data[row][col];
	}

	template<typename T>
	const Vector<T>& Matrix<T>::get_row(const size_t ind) const
	{
		return m_data[ind];
	}

	template<typename T>
	Vector<T> Matrix<T>::get_col(const size_t ind) const
	{
		Vector<T> new_vec = Vector<T>();
		for (size_t i = 0; i < this->shape().first; ++i)
			new_vec.push_back(this->at(i, ind));
		return new_vec;
	}

	template<typename T>
	Matrix<T> Matrix<T>::copy_values() const
	{
		Matrix<T> res = Matrix(this->shape().first, this->shape().second, T());
		for (size_t i = 0; i < this->shape().first; ++i)
			for (size_t j = 0; j < this->shape().second; ++j)
				res.at(i, j) = this->at(i, j);
		return res;
	}

	template<typename T>
	void Matrix<T>::print() const
	{
		for (size_t i = 0; i < this->shape().first; ++i) {
			for (size_t j = 0; j < this->shape().second; ++j)
				std::cout << this->at(i, j) << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<typename T>
	Matrix<T> identity_matrix(size_t size)
	{
		Matrix<T> new_mat(size, size, 0);
		for (size_t i = 0; i < size; ++i)
			new_mat.at(i, i) = T(1);
		return new_mat;
	}
}