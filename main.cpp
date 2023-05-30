#pragma once

#include "Random.h"
#include <iostream>


int main()
{
    for (size_t i = 0; i < 100; i++)
        std::cout << kana::random_uniform() * 100 << std::endl;
    std::string val;
    std::cin >> val;
}