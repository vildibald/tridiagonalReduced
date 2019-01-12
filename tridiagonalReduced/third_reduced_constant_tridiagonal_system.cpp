#include "pch.h"
#include "third_reduced_constant_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* third_reduced_constant_tridiagonal_system::clone_impl() const
{
	return new third_reduced_constant_tridiagonal_system(*this);
}

third_reduced_constant_tridiagonal_system::third_reduced_constant_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(std::move(tridiagonal_system), 8)
{
}

void third_reduced_constant_tridiagonal_system::prepare()
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

	const auto a8 = a4 * a4;
	const auto bbyyy = 2 * a4 * c4 - bbxxx * bbxxx;
	const auto c8 = c4 * c4;

#pragma omp parallel for
	for (int index = 0; index < ec; ++index)
	{
		const auto i = 8 * index + 7;

		ld[index] = a8;

		md[index] = bbyyy;

		ud[index] = c8;

		const auto rho1 = r[i - 7];
		const auto rho2 = r[i - 6];
		const auto rho3 = r[i - 5];
		const auto rho4 = r[i - 4];
		const auto rho5 = r[i - 3];
		const auto rho6 = r[i - 2];
		const auto rho7 = r[i - 1];
		const auto rho8 = r[i];
		const auto rho9 = r[i + 1];
		const auto rho10 = r[i + 2];
		const auto rho11 = r[i + 3];
		const auto rho12 = r[i + 4];
		const auto rho13 = r[i + 5];
		const auto rho14 = r[i + 6];
		const auto rho15 = r[i + 7];

		const auto rabc1 = a * rho1 - b * rho2 + c * rho3;
		const auto rabc3 = a * rho3 - b * rho4 + c * rho5;
		const auto rabc5 = a * rho5 - b * rho6 + c * rho7;
		const auto rabc7 = a * rho7 - b * rho8 + c * rho9;
		const auto rabc9 = a * rho9 - b * rho10 + c * rho11;
		const auto rabc11 = a * rho11 - b * rho12 + c * rho13;
		const auto rabc13 = a * rho13 - b * rho14 + c * rho15;

		const auto raabbcc1 = aa * rabc1 - two_ac_min_bb * rabc3 + cc * rabc5;
		const auto raabbcc5 = aa * rabc5 - two_ac_min_bb * rabc7 + cc * rabc9;
		const auto raabbcc9 = aa * rabc9 - two_ac_min_bb * rabc11 + cc * rabc13;

		rhs[index + 1] = a4 * raabbcc1 -
			bbxxx * raabbcc5 +
			c4 * raabbcc9;
	}

	rhs[0] = df;
	rhs[1] -= a8 * df;
	rhs[rhs.size() - 2] -= c8 * dl;
	rhs[rhs.size() - 1] = dl;
}

std::vector<number> third_reduced_constant_tridiagonal_system::compute_rest()
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

	const auto a4 = aa * aa;
	const auto bbxxx = 2 * aa * cc - two_ac_min_bb * two_ac_min_bb;
	const auto c4 = cc * cc;

	const auto one_div_bbxxx = 1.0 / bbxxx;
	const auto one_div_two_ac_min_bb = 1.0 / two_ac_min_bb;
	const auto one_div_b = 1.0 / b;

#pragma omp parallel for
	for (int index = 0; index < ec + 1; ++index)
	{
		const auto i = 8 * index + 7;

		const auto rho1 = r[i - 7];
		const auto rho2 = r[i - 6];
		const auto rho3 = r[i - 5];
		const auto rho4 = r[i - 4];
		const auto rho5 = r[i - 3];
		const auto rho6 = r[i - 2];
		const auto rho7 = r[i - 1];

		const auto rabc1 = a * rho1 - b * rho2 + c * rho3;
		const auto rabc5 = a * rho5 - b * rho6 + c * rho7;

		const auto d0 = rhs[index];
		const auto d8 = rhs[index + 1];

		const auto d4 = one_div_bbxxx * (
			aa * rabc1 -
			two_ac_min_bb * (a * rho3 - b * rho4 + c * rho5) +
			cc * rabc5 -
			a4 * d0 - c4 * d8);

		const auto d2 = one_div_two_ac_min_bb * (
			rabc1
			- aa * d0 - cc * d4
			);

		const auto d6 = one_div_two_ac_min_bb * (
			rabc5
			- aa * d4 - cc * d8
			);

		result[i - 7] = d0;
		result[i - 6] = (rho1 - a * d0 - c * d2) * one_div_b;
		result[i - 5] = d2;
		result[i - 4] = (rho3 - a * d2 - c * d4) * one_div_b;
		result[i - 3] = d4;
		result[i - 2] = (rho5 - a * d4 - c * d6) * one_div_b;
		result[i - 1] = d6;
		result[i] = (rho7 - a * d6 - c * d8) * one_div_b;
	}
	result.back() = dl;

	timer.stop();

	return result;
}