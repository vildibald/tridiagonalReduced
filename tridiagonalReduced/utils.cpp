#include "pch.h"
#include "utils.h"
#include <omp.h>
#include <mkl.h>

namespace utils
{
	int logical_thread_count = std::thread::hardware_concurrency();

	//	void solve_tridiagonal_system(number* lower_diagonal,
	//	                              number* main_diagonal,
	//	                              number* upper_diagonal,
	//	                              number* right_side,
	//	                              const int num_equations)
	//	{
	//		const auto buffer = new number[num_equations];
	//		solve_tridiagonal_system_buffered(
	//			lower_diagonal,
	//			main_diagonal,
	//			upper_diagonal,
	//			right_side,
	//			num_equations,
	//			buffer);
	//		delete[] buffer;
	//	}
	//
	//	void solve_tridiagonal_system_buffered(const number* lower_diagonal,
	//	                                       const number* main_diagonal,
	//	                                       number* upper_diagonal,
	//	                                       number* right_side,
	//	                                       const int num_equations,
	//	                                       number* buffer)
	//	{
	//		std::copy(upper_diagonal, upper_diagonal + num_equations, buffer);
	//		//		memcpy(buffer, upperDiagonal, numEquations);
	//		buffer[0] /= main_diagonal[0];
	//		right_side[0] /= main_diagonal[0];
	//		for (int i = 1; i < num_equations; i++)
	//		{
	//			const auto m = 1 / (main_diagonal[i] - lower_diagonal[i] *
	//				buffer[i - 1]);
	//			buffer[i] *= m;
	//			right_side[i] = (right_side[i] - lower_diagonal[i] *
	//				right_side[i - 1]) * m;
	//		}
	//		for (int i = num_equations - 1; i-- > 0;)
	//		{
	//			right_side[i] -= buffer[i] * right_side[i + 1];
	//		}
	//	}

	void solve_deboor_tridiagonal_system(const number main_diagonal_value,
		number* right_side,
		const int num_equations,
		const number last_main_diagonal_value)
	{
		const auto buffer = new number[num_equations];
		solve_deboor_tridiagonal_system_buffered(main_diagonal_value, right_side,
			num_equations, buffer, last_main_diagonal_value);
		delete[] buffer;
	}

	void solve_deboor_tridiagonal_system(const number lower_diagonal_value,
		const number main_diagonal_value,
		const number upper_diagonal_value,
		number* right_side,
		const int num_equations,
		const number last_main_diagonal_value)
	{
		const auto buffer = new number[num_equations];
		solve_deboor_tridiagonal_system_buffered(lower_diagonal_value,
			main_diagonal_value,
			upper_diagonal_value,
			right_side,
			num_equations,
			buffer,
			last_main_diagonal_value);
		delete[] buffer;
	}

	void solve_deboor_tridiagonal_system_buffered(const number lower_diagonal_value,
		const number main_diagonal_value,
		const number upper_diagonal_value,
		number* right_side,
		const int num_equations,
		number* buffer,
		number last_main_diagonal_value)
	{
		if (last_main_diagonal_value == DBL_MIN)
			last_main_diagonal_value = main_diagonal_value;
		auto m0 = 1 / main_diagonal_value;
		buffer[0] = upper_diagonal_value * m0;
		right_side[0] *= m0;
		const auto lastindex = num_equations - 1;
		for (int i = 1; i < lastindex; i++)
		{
			const auto m = 1 / (main_diagonal_value - lower_diagonal_value *
				buffer[i - 1]);
			buffer[i] = upper_diagonal_value * m;
			right_side[i] = (right_side[i] - lower_diagonal_value *
				right_side[i - 1]) * m;
		}

		m0 = 1 / (last_main_diagonal_value - lower_diagonal_value *
			buffer[lastindex - 1]);
		buffer[lastindex] = upper_diagonal_value * m0;
		right_side[lastindex] = (right_side[lastindex] - lower_diagonal_value *
			right_side[lastindex - 1]) * m0;

		for (int i = num_equations - 1; i-- > 0;)
		{
			right_side[i] -= buffer[i] * right_side[i + 1];
		}
	}

	void solve_deboor_tridiagonal_system_buffered(const number main_diagonal_value,
		number* right_side,
		const int num_equations,
		number* buffer,
		number last_main_diagonal_value)
	{
		if (last_main_diagonal_value == DBL_MIN)
			last_main_diagonal_value = main_diagonal_value;
		auto m0 = 1 / main_diagonal_value;
		buffer[0] = m0;
		right_side[0] *= m0;
		const auto lastindex = num_equations - 1;
		for (int i = 1; i < lastindex; i++)
		{
			auto m = 1 / (main_diagonal_value - buffer[i - 1]);
			buffer[i] = m;
			right_side[i] = (right_side[i] - right_side[i - 1]) * m;
		}

		m0 = 1 / (last_main_diagonal_value - buffer[lastindex - 1]);
		buffer[lastindex] = m0;
		right_side[lastindex] = (right_side[lastindex]
			- right_side[lastindex - 1]) * m0;

		for (int i = num_equations - 1; i-- > 0;)
		{
			right_side[i] -= buffer[i] * right_side[i + 1];
		}
	}

	void solve_tridiagonal_system_buffered_mkl(number* lower_diagonal,
		number* main_diagonal,
		number* upper_diagonal,
		number* right_sides,
		const int num_equations,
		const int num_rhs)
	{
		int factor_info = 0;
		sdttrfb(&num_equations,
			lower_diagonal,
			main_diagonal,
			upper_diagonal,
			&factor_info);

		const char trans = 'N';
		int solve_info = 0;

		sdttrsb(&trans,
			&num_equations,
			&num_rhs,
			lower_diagonal,
			main_diagonal,
			upper_diagonal,
			right_sides,
			&num_equations,
			&solve_info);
	}
}