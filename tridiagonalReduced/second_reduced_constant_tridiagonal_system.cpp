#include "pch.h"
#include "second_reduced_constant_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* second_reduced_constant_tridiagonal_system::clone_impl() const
{
	return new second_reduced_constant_tridiagonal_system(*this);
}

second_reduced_constant_tridiagonal_system::second_reduced_constant_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(std::move(tridiagonal_system), 4)
{
}

void second_reduced_constant_tridiagonal_system::prepare()
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
	const auto bb = b * b;
	const auto two_ac_min_bb = 2 * a * c - bb;
	const auto cc = c * c;

	const auto a4 = aa * aa;
	const auto bbxxx = 2 * aa * cc - two_ac_min_bb * two_ac_min_bb;
	const auto c4 = cc * cc;

#pragma omp parallel for
	for (int index = 0; index < ec; ++index)
	{
		const auto i = 4 * index + 3;

		ld[index] = a4;

		md[index] = bbxxx;

		ud[index] = c4;

		const auto rho1 = r[i - 3];
		const auto rho2 = r[i - 2];
		const auto rho3 = r[i - 1];
		const auto rho4 = r[i];
		const auto rho5 = r[i + 1];
		const auto rho6 = r[i + 2];
		const auto rho7 = r[i + 3];

		const auto rabc1 = a * rho1 - b * rho2 + c * rho3;
		const auto rabc3 = a * rho3 - b * rho4 + c * rho5;
		const auto rabc5 = a * rho5 - b * rho6 + c * rho7;

		rhs[index + 1] = aa * rabc1 - two_ac_min_bb * rabc3 + cc * rabc5;
	}

	rhs[0] = df;
	rhs[1] -= a4 * df;
	rhs[ec] -= c4 * dl;
	rhs[ec + 1] = dl;
}

std::vector<number> second_reduced_constant_tridiagonal_system::compute_rest()
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

	const auto aa = a * a;
	const auto bb = b * b;
	const auto two_ac_min_bb = 2 * a * c - bb;
	const auto cc = c * c;

	const auto one_div_two_ac_min_bb = 1.0 / two_ac_min_bb;
	const auto one_div_b = 1.0 / b;

#pragma omp parallel for
	for (int index = 0; index < ec + 1; ++index)
	{
		const auto i = 4 * index + 3;

		const auto rho1 = r[i - 3];
		const auto rho2 = r[i - 2];
		const auto rho3 = r[i - 1];

		const auto d0 = rhs[index];;
		const auto d4 = rhs[index + 1];

		const auto d2 = one_div_two_ac_min_bb * (
			a * rho1 - b * rho2 + c * rho3 - aa * d0 - cc * d4
			);

		result[i - 3] = d0;
		result[i - 2] = (rho1 - a * d0 - c * d2) * one_div_b;
		result[i - 1] = d2;
		result[i] = (rho3 - a * d2 - c * d4) * one_div_b;
	}
	result.back() = dl;

	timer.stop();

	return result;
}
