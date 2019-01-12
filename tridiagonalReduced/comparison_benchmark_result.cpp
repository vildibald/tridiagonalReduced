#include "pch.h"
#include "comparison_benchmark_result.h"

comparison_benchmark_result& comparison_benchmark_result::add(unsigned long long time)
{
	algorithm_times_.emplace_back(time);
	return *this;
}

unsigned long long comparison_benchmark_result::operator[](const int i) const
{
	return algorithm_times_[i];
}

number comparison_benchmark_result::ratio(const int i, const int j) const
{
	return static_cast<number>(algorithm_times_[i]) / static_cast<number>(algorithm_times_[j]);
}
