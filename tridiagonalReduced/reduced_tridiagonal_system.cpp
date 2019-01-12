#include "pch.h"
#include "reduced_tridiagonal_system.h"
#include "utils.h"

standard_tridiagonal_system& reduced_tridiagonal_system::full_system()
{
	return full_system_;
}

reduced_tridiagonal_system::reduced_tridiagonal_system(standard_tridiagonal_system full_system, int factor):
	tridiagonal_system(std::vector<number>(full_system.equation_count() / factor),
	                   std::vector<number>(full_system.equation_count() / factor),
	                   std::vector<number>(full_system.equation_count() / factor),
	                   std::vector<number>(2 + full_system.equation_count() / factor),
	                   full_system.d_first(),
	                   full_system.d_last()
	),
	full_system_(std::move(full_system))
{
}


reduced_tridiagonal_system::~reduced_tridiagonal_system() = default;

int reduced_tridiagonal_system::reduction_factor() const
{
	return full_system_.equation_count() / equation_count();
}

std::vector<number> reduced_tridiagonal_system::solve()
{
	auto& timer = get_timer();

	prepare();

	timer.start_all_time();

	const auto ec = equation_count();
	auto& ld = lower_diagonal();
	auto& md = main_diagonal();
	auto& ud = upper_diagonal();
	auto& rhs = right_side();

	timer.start_execution_time();

	utils::solve_tridiagonal_system_buffered_mkl(
		&ld[1],
		&md[0],
		&ud[0],
		&rhs[1],
		ec,
		1
	);

	timer.stop();

	return compute_rest();
}
