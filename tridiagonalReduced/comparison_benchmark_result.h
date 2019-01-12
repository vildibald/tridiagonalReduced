#pragma once
#include <vector>
#include "number.h"

class comparison_benchmark_result {
	std::vector<unsigned long long> algorithm_times_;

public:
	comparison_benchmark_result& add(unsigned long long time);

	unsigned long long operator[](int i) const;

	number ratio(int i, int j) const;
};
