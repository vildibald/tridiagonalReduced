#include "pch.h"
#include "first_reduced_constant_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* first_reduced_constant_tridiagonal_system::clone_impl() const
{
	return new first_reduced_constant_tridiagonal_system(*this);
}

first_reduced_constant_tridiagonal_system::first_reduced_constant_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(std::move(tridiagonal_system), 2)
{
}

void first_reduced_constant_tridiagonal_system::prepare()
{
	auto& full = full_system();
	auto& timer = get_timer();

	timer.start();

	const auto a = full.lower_diagonal()[0];
	const auto b = full.main_diagonal()[0];
	const auto c = full.upper_diagonal()[0];
	const auto& r = full.right_side();
	const auto ec = equation_count();
	auto& ld = lower_diagonal();
	auto& md = main_diagonal();
	auto& ud = upper_diagonal();
	auto& rhs = right_side();
	const auto df = d_first();
	const auto dl = d_last();

	const auto aa = a * a;
	const auto two_ac_min_bb = 2 * a * c - b * b;
	const auto cc = c * c;

#pragma omp parallel for
	for (int index = 0; index < ec; ++index)
	{
		const auto i = 2 * index + 1;

		ld[index] = aa;

		md[index] =
			two_ac_min_bb;

		ud[index] = cc;

		rhs[index + 1] =
			a * r[i - 1] - b * r[i] + c * r[i + 1];
	}
	rhs[0] = df;
	rhs[1] -= aa * df;
	rhs[ec] -= cc * dl;
	rhs[ec + 1] = dl;
}

std::vector<number> first_reduced_constant_tridiagonal_system::compute_rest()
{
	auto& full = full_system();
	auto& timer = get_timer();

	timer.start_all_time();

	std::vector<number> result(full.equation_count() + 2);

	timer.start_execution_time();

	const auto a = full.lower_diagonal()[0];
	const auto b = full.main_diagonal()[0];
	const auto c = full.upper_diagonal()[0];
	const auto& r = full.right_side();
	const auto ec = equation_count();
	auto& rhs = right_side();
	const auto dl = d_last();

	const auto one_div_md = 1.0 / b;

#pragma omp parallel for
	for (int index = 0; index < ec + 1; ++index)
	{
		const auto i = 2 * index + 1;

		const auto dp = rhs[index];
		const auto dn = rhs[index + 1];

		result[i - 1] = dp;
		result[i] = (r[i - 1] - a * dp - c * dn) * one_div_md;
	}
	result.back() = dl;

	timer.stop();

	return result;
}