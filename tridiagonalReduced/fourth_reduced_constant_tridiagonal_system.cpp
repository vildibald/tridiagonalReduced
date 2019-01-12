#include "pch.h"
#include "fourth_reduced_constant_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* fourth_reduced_constant_tridiagonal_system::clone_impl() const
{
	return new fourth_reduced_constant_tridiagonal_system(*this);
}

fourth_reduced_constant_tridiagonal_system::fourth_reduced_constant_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(
		std::move(tridiagonal_system), 16)
{
}

void fourth_reduced_constant_tridiagonal_system::prepare()
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

	const auto a16 = a8 * a8;
	const auto bbzzz = 2 * a8 * c8 - bbyyy * bbyyy;
	const auto c16 = c8 * c8;

#pragma omp parallel for
	for (int index = 0; index < ec; ++index)
	{
		const auto i = 16 * index + 15;

		ld[index] = a16;

		md[index] = bbzzz;

		ud[index] = c16;

		const auto rho1 = r[i - 15];
		const auto rho2 = r[i - 14];
		const auto rho3 = r[i - 13];
		const auto rho4 = r[i - 12];
		const auto rho5 = r[i - 11];
		const auto rho6 = r[i - 10];
		const auto rho7 = r[i - 9];
		const auto rho8 = r[i - 8];
		const auto rho9 = r[i - 7];
		const auto rho10 = r[i - 6];
		const auto rho11 = r[i - 5];
		const auto rho12 = r[i - 4];
		const auto rho13 = r[i - 3];
		const auto rho14 = r[i - 2];
		const auto rho15 = r[i - 1];
		const auto rho16 = r[i];
		const auto rho17 = r[i + 1];
		const auto rho18 = r[i + 2];
		const auto rho19 = r[i + 3];
		const auto rho20 = r[i + 4];
		const auto rho21 = r[i + 5];
		const auto rho22 = r[i + 6];
		const auto rho23 = r[i + 7];
		const auto rho24 = r[i + 8];
		const auto rho25 = r[i + 9];
		const auto rho26 = r[i + 10];
		const auto rho27 = r[i + 11];
		const auto rho28 = r[i + 12];
		const auto rho29 = r[i + 13];
		const auto rho30 = r[i + 14];
		const auto rho31 = r[i + 15];

		const auto rabc1 = a * rho1 - b * rho2 + c * rho3;
		const auto rabc3 = a * rho3 - b * rho4 + c * rho5;
		const auto rabc5 = a * rho5 - b * rho6 + c * rho7;
		const auto rabc7 = a * rho7 - b * rho8 + c * rho9;
		const auto rabc9 = a * rho9 - b * rho10 + c * rho11;
		const auto rabc11 = a * rho11 - b * rho12 + c * rho13;
		const auto rabc13 = a * rho13 - b * rho14 + c * rho15;
		const auto rabc15 = a * rho15 - b * rho16 + c * rho17;
		const auto rabc17 = a * rho17 - b * rho18 + c * rho19;
		const auto rabc19 = a * rho19 - b * rho20 + c * rho21;
		const auto rabc21 = a * rho21 - b * rho22 + c * rho23;
		const auto rabc23 = a * rho23 - b * rho24 + c * rho25;
		const auto rabc25 = a * rho25 - b * rho26 + c * rho27;
		const auto rabc27 = a * rho27 - b * rho28 + c * rho29;
		const auto rabc29 = a * rho29 - b * rho30 + c * rho31;

		const auto raabbcc1 = aa * rabc1 - two_ac_min_bb * rabc3 + cc * rabc5;
		const auto raabbcc5 = aa * rabc5 - two_ac_min_bb * rabc7 + cc * rabc9;
		const auto raabbcc9 = aa * rabc9 - two_ac_min_bb * rabc11 + cc * rabc13;
		const auto raabbcc13 = aa * rabc13 - two_ac_min_bb * rabc15 + cc * rabc17;
		const auto raabbcc17 = aa * rabc17 - two_ac_min_bb * rabc19 + cc * rabc21;
		const auto raabbcc21 = aa * rabc21 - two_ac_min_bb * rabc23 + cc * rabc25;
		const auto raabbcc25 = aa * rabc25 - two_ac_min_bb * rabc27 + cc * rabc29;

		const auto ra4b4c41 = a4 * raabbcc1 - bbxxx * raabbcc5 + c4 * raabbcc9;
		const auto ra4b4c49 = a4 * raabbcc9 - bbxxx * raabbcc13 + c4 * raabbcc17;
		const auto ra4b4c417 = a4 * raabbcc17 - bbxxx * raabbcc21 + c4 * raabbcc25;


		rhs[index + 1] = a8 * ra4b4c41 -
			bbyyy * ra4b4c49 +
			c8 * ra4b4c417;
	}

	rhs[0] = df;
	rhs[1] -= a16 * df;
	rhs[ec] -= c16 * dl;
	rhs[ec + 1] = dl;

	utils::solve_tridiagonal_system_buffered_mkl(
		&ld[1],
		&md[0],
		&ud[0],
		&rhs[1],
		ec,
		1
	);
}

std::vector<number> fourth_reduced_constant_tridiagonal_system::compute_rest()
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

	const auto bbyyy = 2 * a4 * c4 - bbxxx * bbxxx;

	const auto one_div_bbyyy = 1.0 / bbyyy;
	const auto one_div_bbxxx = 1.0 / bbxxx;
	const auto one_div_two_ac_min_bb = 1.0 / two_ac_min_bb;
	const auto one_div_b = 1.0 / b;

#pragma omp parallel for
	for (int index = 0; index < ec + 1; ++index)
	{
		const auto i = 16 * index + 15;

		const auto rho1 = r[i - 15];
		const auto rho2 = r[i - 14];
		const auto rho3 = r[i - 13];
		const auto rho4 = r[i - 12];
		const auto rho5 = r[i - 11];
		const auto rho6 = r[i - 10];
		const auto rho7 = r[i - 9];
		const auto rho8 = r[i - 8];
		const auto rho9 = r[i - 7];
		const auto rho10 = r[i - 6];
		const auto rho11 = r[i - 5];
		const auto rho12 = r[i - 4];
		const auto rho13 = r[i - 3];
		const auto rho14 = r[i - 2];
		const auto rho15 = r[i - 1];

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

		const auto d0 = rhs[index];
		const auto d16 = rhs[index + 1];


		const auto d8 = one_div_bbyyy * (a4 * raabbcc1 -
			bbxxx * raabbcc5 +
			c4 * raabbcc9 - d0 - d16);

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

		const auto d12 = one_div_bbxxx * (
			aa * rabc9 -
			two_ac_min_bb * rabc11 +
			cc * rabc13 -
			d8 - d16);

		const auto d10 = one_div_two_ac_min_bb * (
			rabc9
			- d8 - d12
			);

		const auto d14 = one_div_two_ac_min_bb * (
			rabc13
			- d12 - d16
			);

		result[i - 15] = d0;
		result[i - 14] = (rho1 - a * d0 - c * d2) * one_div_b;
		result[i - 13] = d2;
		result[i - 12] = (rho3 - a * d2 - c * d4) * one_div_b;
		result[i - 11] = d4;
		result[i - 10] = (rho5 - a * d4 - c * d6) * one_div_b;
		result[i - 9] = d6;
		result[i - 8] = (rho7 - a * d6 - c * d8) * one_div_b;
		result[i - 7] = d8;
		result[i - 6] = (rho9 - a * d8 - c * d10) * one_div_b;
		result[i - 5] = d10;
		result[i - 4] = (rho11 - a * d10 - c * d12) * one_div_b;
		result[i - 3] = d12;
		result[i - 2] = (rho13 - a * d12 - c * d14) * one_div_b;
		result[i - 1] = d14;
		result[i] = (rho15 - a * d14 - c * d16) * one_div_b;
	}
	result.back() = dl;

	timer.stop();

	return result;
}