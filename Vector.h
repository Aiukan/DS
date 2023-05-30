#pragma once

#include <vector>
#include <functional>

namespace kana
{
	template<typename T>
	class Vector
	{
	public:
		Vector();
		Vector(const std::vector<T> vec);
		Vector(const std::initializer_list<T> list);
		Vector(const size_t size, const T val);
		Vector(const size_t size, const std::function<T()> f);
		virtual ~Vector() = default;

		Vector<T>(const Vector<T>& other) = default;
		Vector<T>& operator=(const Vector<T>& other) = default;

		Vector<T>(Vector<T>&& other) noexcept;
		Vector<T>& operator=(Vector<T>&& other) noexcept;

		Vector<T>& operator-();

		Vector<T> operator+(const Vector<T>& other) const;
		Vector<T>& operator+=(const Vector<T>& other);

		Vector<T> operator-(const Vector<T>& other) const;
		Vector<T>& operator-=(const Vector<T>& other);

		Vector<T> operator*(const T val) const;
		Vector<T>& operator*=(const T val);

		Vector<T> operator/(const T val) const;
		Vector<T>& operator/=(const T val);

		T& operator[](const size_t ind);

		T operator[](const size_t ind) const;

		Vector<T> copy_values() const;
		void sort();
		Vector<T> sorted() const;
		void push_back(const T val);
		void insert(const size_t index, const T val);
		T remove(const size_t index);


		size_t size() const;

		void print() const;
	private:
		std::vector<T> m_data;
	};

	template<typename T>
	std::_Vector_iterator<T> begin(const Vector<T> vec);

	template<typename T>
	std::_Vector_iterator<T> end(const Vector<T> vec);

	template<typename T>
	std::_Vector_const_iterator<T> cbegin(const Vector<T> vec);

	template<typename T>
	std::_Vector_const_iterator<T> cend(const Vector<T> vec);

}