#pragma once

#include "Vector.h"

namespace kana
{
	double random_uniform(const double left = 0, const double right = 1);

	double random_normal(const double mu = 0, const double sigma = 1);

	template<typename T>
	void shuffle(Vector<T>& vec);

	template<typename T>
	Vector<T> shuffled(const Vector<T>& vec);
}