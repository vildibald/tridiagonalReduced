#pragma once
#include <windows.h>
#include <chrono>

class stop_watch
{
	std::chrono::steady_clock::time_point start_{};
	std::chrono::steady_clock::time_point stop_{};
	long long elapsed_ = 0;

public:
	void start()
	{
		start_ = std::chrono::high_resolution_clock::now();
	}

	void stop()
	{
		stop_ = std::chrono::high_resolution_clock::now();
		const auto diff = stop_ - start_;
		elapsed_ += std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
	}

	long long elapsed_time() const
	{
		return elapsed_;
	}

	void reset()
	{
		elapsed_ = 0;
	}

	stop_watch& operator+=(stop_watch& other)
	{
		elapsed_ += other.elapsed_;
		return *this;
	}
};