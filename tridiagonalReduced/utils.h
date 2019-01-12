#pragma once

#include <omp.h>
#include <thread>
#include <vector>
#include "number.h"

namespace utils
{
	extern int logical_thread_count;

	template <typename Iterator, typename Function>
	void parallelizable_for(Iterator from, Iterator to, Iterator increment_by, const bool in_parallel,
		Function& function)
	{
		if (in_parallel)
		{
#pragma omp parallel for
			for (Iterator i = from; i < to; i += increment_by)
			{
				function(i);
			}
		}
		else
		{
			// Loop for sequential computations.
			// '#pragma omp parallel for if(in_parallel)' in case of in_parallel==false will execute
			// such loop in one thread, but still with overhead of OpenMP thread creation.
			for (Iterator i = from; i < to; i += increment_by)
			{
				function(i);
			}
		}
	}

	//	void solve_tridiagonal_system(number* lower_diagonal,
	//	                              number* main_diagonal,
	//	                              number* upper_diagonal,
	//	                              number* right_side,
	//	                              int num_equations);

	//	void solve_tridiagonal_system_buffered(number* lower_diagonal,
	//	                                       number* main_diagonal,
	//	                                       number* upper_diagonal,
	//	                                       number* right_side,
	//	                                       int num_equations,
	//	                                       number* buffer);

	void solve_deboor_tridiagonal_system(number main_diagonal_value,
		number* right_side,
		int num_equations,
		number last_main_diagonal_value = DBL_MIN);

	void solve_deboor_tridiagonal_system(number lower_diagonal_value,
		number main_diagonal_value,
		number upper_diagonal_value,
		number* right_side,
		int num_equations,
		number last_main_diagonal_value = DBL_MIN);

	void solve_deboor_tridiagonal_system_buffered(number lower_diagonal_value,
		number main_diagonal_value,
		number upper_diagonal_value,
		number* right_side,
		int num_equations,
		number* buffer,
		number last_main_diagonal_value = DBL_MIN);

	void solve_deboor_tridiagonal_system_buffered(number main_diagonal_value,
		number* right_side,
		int num_equations,
		number* buffer,
		number last_main_diagonal_value = DBL_MIN);

	void solve_tridiagonal_system_buffered_mkl(number* lower_diagonal,
		number* main_diagonal,
		number* upper_diagonal,
		number* right_sides,
		int num_equations,
		int num_rhs);
};

