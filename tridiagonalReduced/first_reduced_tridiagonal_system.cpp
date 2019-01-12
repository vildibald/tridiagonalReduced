#include "pch.h"
#include "first_reduced_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* first_reduced_tridiagonal_system::clone_impl() const
{
	return new first_reduced_tridiagonal_system(*this);
}

first_reduced_tridiagonal_system::first_reduced_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(
	std::move(tridiagonal_system), 2)
{
}

void first_reduced_tridiagonal_system::prepare()
{
	auto& full = full_system();
	auto& timer = get_timer();

	timer.start();

	const auto& a = full.lower_diagonal();
	const auto& b = full.main_diagonal();
	const auto& c = full.upper_diagonal();
	const auto& r = full.right_side();
	const auto ec = equation_count();
	auto& ld = lower_diagonal();
	auto& md = main_diagonal();
	auto& ud = upper_diagonal();
	auto& rhs = right_side();
	const auto df = d_first();
	const auto dl = d_last();

#pragma omp parallel for
	for (int index = 0; index < ec; ++index)
	{
		const auto i = 2 * index + 1;

		const auto a1 = a[i - 1];
		const auto a2 = a[i];
		const auto a3 = a[i + 1];

		const auto b1 = b[i - 1];
		const auto b2 = b[i];
		const auto b3 = b[i + 1];

		const auto c1 = c[i - 1];
		const auto c2 = c[i];
		const auto c3 = c[i + 1];

		const auto rho1 = r[i - 1];
		const auto rho2 = r[i];
		const auto rho3 = r[i + 1];

		ld[index] = b3 * a2 * a1;

		md[index] =
			b3 * (a2 * c1 - b1 * b2) + b1 * c2 * a3;

		ud[index] = b1 * c2 * c3;

		rhs[index + 1] =
			b3 * (a2 * rho1 - b1 * rho2) + b1 * c2 * rho3;
	}
	rhs[0] = df;
	rhs[1] -= ld[0] * df;
	rhs[ec] -= ud.back() * dl;
	rhs[ec + 1] = dl;

	timer.stop();
}

std::vector<number> first_reduced_tridiagonal_system::compute_rest()
{
	auto& timer = get_timer();
	timer.start_all_time();

	auto& full = full_system();
	
	std::vector<number> result(full.equation_count() + 2);

	timer.start_execution_time();

	const auto& a = full.lower_diagonal();
	const auto& b = full.main_diagonal();
	const auto& c = full.upper_diagonal();
	const auto& r = full.right_side();
	const auto ec = equation_count();
	auto& rhs = right_side();
	const auto dl = d_last();

#pragma omp parallel for
	for (int index = 0; index < ec + 1; ++index)
	{
		const auto i = 2 * index + 1;

		const auto dp = rhs[index];
		const auto dn = rhs[index + 1];

		result[i - 1] = dp;
		result[i] =
			(r[i - 1] - a[i - 1] * dp - c[i - 1] * dn) * (1 / b[i - 1]);
	}
	result.back() = dl;

	timer.stop();

	return result;
}