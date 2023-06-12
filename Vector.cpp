#pragma once

#include "Vector.h"
#include "Random.h"
#include <vector>
#include <iostream>
#include <algorithm>

namespace kana 
{
	template<typename T>
	Vector<T>::Vector()
		: m_data()
	{
	}

	template<typename T>
	Vector<T>::Vector(const std::vector<T> vec)
		: m_data(vec)
	{
	}

	template<typename T>
	Vector<T>::Vector(const std::initializer_list<T> list)
		: m_data(list)
	{
	}

	template<typename T>
	Vector<T>::Vector(const size_t size, const T val)
		: m_data(size, val)
	{
	}

	template<typename T>
	Vector<T>::Vector(const size_t size, const std::function<T()> f)
		: m_data(size)
	{
		std::generate(begin(m_data), end(m_data), f);
	}

	template<typename T>
	Vector<T>::Vector(Vector<T>&& other) noexcept
		: m_data(std::move(other.m_data))
	{
	}

	template<typename T>
	Vector<T>& Vector<T>::operator=(Vector<T>&& other) noexcept
	{
		this->m_data = move(other.m_data);
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator-()
	{
		for (size_t i = 0; i < this->size(); ++i)
		{
			m_data[i] = -m_data[i];
		}

		return *this;
	}

	template<typename T>
	Vector<T> Vector<T>::operator+(const Vector<T>& other) const
	{
		Vector<T> new_vec(this->size(), 0);
		for (size_t i = 0; i < this->size(); ++i)
		{
			new_vec[i] = this->m_data[i] + other.m_data[i];
		}
		return new_vec;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator+=(const Vector<T>& other)
	{
		for (size_t i = 0; i < other.size(); ++i)
		{
			this->m_data[i] += other.m_data[i];
		}
		return *this;
	}

	template<typename T>
	Vector<T> Vector<T>::operator-(const Vector<T>& other) const
	{
		Vector<T> new_vec(this->size(), 0);
		for (size_t i = 0; i < this->size(); ++i)
		{
			new_vec[i] = this->m_data[i] - other.m_data[i];
		}
		return new_vec;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator-=(const Vector<T>& other)
	{
		for (size_t i = 0; i < other.size(); ++i)
		{
			this->m_data[i] -= other.m_data[i];
		}
		return *this;
	}

	template<typename T>
	Vector<T> Vector<T>::operator*(const T val) const
	{
		Vector<T> new_vec(this->size(), val);
		for (size_t i = 0; i < this->size(); ++i)
		{
			new_vec[i] *= m_data[i];
		}
		return new_vec;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator*=(const T val)
	{
		for (size_t i = 0; i < this->size(); ++i)
		{
			m_data[i] *= val;
		}
		return *this;
	}

	template<typename T>
	Vector<T> Vector<T>::operator/(const T val) const
	{
		Vector<T> new_vec(this->size(), 0);
		for (size_t i = 0; i < this->size(); ++i)
		{
			new_vec[i] = m_data[i] / val;
		}
		return new_vec;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator/=(const T val)
	{
		for (size_t i = 0; i < this->size(); ++i)
		{
			m_data[i] /= val;
		}
		return *this;
	}

	template<typename T>
	T& Vector<T>::operator[](const size_t ind)
	{
		return m_data[ind];
	}

	template<typename T>
	T Vector<T>::operator[](const size_t ind) const
	{
		return m_data[ind];
	}

	template<typename T>
	size_t Vector<T>::size() const
	{
		return m_data.size();
	}

	template<typename T>
	void Vector<T>::sort()
	{
		std::sort(begin(m_data), end(m_data));
	}

	template<typename T>
	Vector<T> Vector<T>::sorted() const
	{
		Vector<T> new_vec(*this);
		new_vec.sort();
		return new_vec;
	}

	template<typename T>
	void Vector<T>::print() const
	{
		for (size_t i = 0; i < this->size(); ++i)
			std::cout << (*this)[i] << ' ';
		std::cout << std::endl;
	}

	template<typename T>
	void Vector<T>::push_back(const T value)
	{
		m_data.push_back(value);
	}

	template<typename T>
	void Vector<T>::insert(const size_t index, const T value)
	{
		m_data.insert(m_data.begin() + index, value);
	}

	template<typename T>
	T Vector<T>::remove(const size_t index)
	{
		T value = (*this)[index];
		m_data.erase(m_data.begin() + index);
		return value;
	}

	template<typename T>
	Vector<T> Vector<T>::copy_values() const
	{
		Vector<T> target = Vector();
		for (size_t i = 0; i < this->size(); ++i)
			target.push_back((*this)[i]);
		return target;

	}

	template<typename T>
	std::_Vector_iterator<T> begin(const Vector<T> vec)
	{
		return std::begin(vec);
	}

	template<typename T>
	std::_Vector_iterator<T> end(const Vector<T> vec)
	{
		return std::end(vec);
	}

	template<typename T>
	std::_Vector_const_iterator<T> cbegin(const Vector<T> vec)
	{
		return std::cbegin(vec);
	}

	template<typename T>
	std::_Vector_const_iterator<T> cend(const Vector<T> vec)
	{
		return std::cend(vec);
	}
}