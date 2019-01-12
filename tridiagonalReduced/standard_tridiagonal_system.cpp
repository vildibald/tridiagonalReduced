#include "pch.h"
#include "standard_tridiagonal_system.h"
#include <string>
#include "utils.h"

standard_tridiagonal_system::standard_tridiagonal_system(std::vector<number> lower_diagonal,
	std::vector<number> main_diagonal, std::vector<number> upper_diagonal, std::vector<number> right_side,
	number d_first, number d_last, timer timer) : tridiagonal_system(
		std::move(lower_diagonal),
		std::move(main_diagonal), 
		std::move(upper_diagonal), 
		std::move(right_side), 
		d_first, 
		d_last, 
		std::move(timer))
{
}

standard_tridiagonal_system::standard_tridiagonal_system(std::vector<number> lower_diagonal,
	std::vector<number> main_diagonal, std::vector<number> upper_diagonal, std::vector<number> right_side,
	number d_first, number d_last) : tridiagonal_system(
		std::move(lower_diagonal),
		std::move(main_diagonal),
		std::move(upper_diagonal),
		std::move(right_side),
		d_first,
		d_last)
{
}


std::vector<number> standard_tridiagonal_system::solve()
{
	auto& timer = get_timer();
	timer.start_all_time();

	const int ec = equation_count();
	auto& ld = lower_diagonal();
	auto& md = main_diagonal();
	auto& ud = upper_diagonal();
	auto& rhs = right_side();
	const auto df = d_first();
	const auto dl = d_last();

	std::vector<number> result(rhs.size() + 2);

	timer.start_execution_time();

	rhs[0] -= ld[0] * df;
	rhs.back() -= ud.back() * dl;
	utils::solve_tridiagonal_system_buffered_mkl(
		&ld[1],
		&md.front(),
		&ud.front(),
		&rhs.front(),
		ec,
		1
	);

	result[0] = df;
	result.back() = dl;

	std::copy(rhs.begin(), rhs.end(), ++result.begin());

	timer.stop();
	return result;
}
