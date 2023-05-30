#pragma once

#include "Random.h"
#include <random>

namespace kana
{
	std::mt19937 _eng(static_cast<std::mt19937::result_type>(std::time(nullptr)));


	double random_uniform(const double left, const double right)
	{
		static std::uniform_real_distribution<double> dist_uniform;
		return left + dist_uniform(_eng) * (right - left); 

	}

	double random_normal(const double mu, const double sigma)
	{
		static std::normal_distribution<> dist_normal{ 0, 1 };
		return dist_normal(_eng) * sigma + mu;
	}

	template<typename T>
	void shuffle(Vector<T>& vec)
	{
		size_t len = vec.size();
		Vector<T> res = Vector();
		size_t index;
		while (len != 0) {
			index = static_cast<size_t>(random_uniform() * len);
			res.push_back(vec.remove(index));
			len--;
		}
		(*this) = res;
	}

	template<typename T>
	Vector<T> shuffled(const Vector<T>& vec)
	{
		size_t len = this->size();
		Vector<T> current = this->copy_values();
		Vector<T> res = Vector();
		size_t index;
		while (len != 0) {
			index = static_cast<size_t>(random_uniform() * len);
			res.push_back(current.remove(index));
			len--;
		}
		return res;
	}

}