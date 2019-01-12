#pragma once
#include "StopWatch.h"

class timer final {
	stop_watch execution_watch_;
	stop_watch all_watch_;
public:
	long long execution_time() const
	{
		return execution_watch_.elapsed_time();
	}

	long long all_time() const
	{
		return all_watch_.elapsed_time();
	}

	void start_execution_time()
	{
		execution_watch_.start();
	}

	void start_all_time()
	{
		all_watch_.start();
	}

	void stop_execution_time()
	{
		execution_watch_.stop();
	}

	void stop_all_time()
	{
		all_watch_.stop();
	}

	void start()
	{
		all_watch_.start();
		execution_watch_.start();
	}

	void stop()
	{
		execution_watch_.stop();
		all_watch_.stop();
	}

	void reset()
	{
		execution_watch_.reset();
		all_watch_.reset();
	}

	timer& operator+=(timer& other){
		all_watch_ += other.all_watch_;
		execution_watch_ += other.execution_watch_;
		return *this;
	}
};

