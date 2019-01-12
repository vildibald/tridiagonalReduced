#include "pch.h"
#include "second_reduced_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* second_reduced_tridiagonal_system::clone_impl() const
{
	return new second_reduced_tridiagonal_system(*this);
}

second_reduced_tridiagonal_system::second_reduced_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(std::move(tridiagonal_system), 4)
{
}

void second_reduced_tridiagonal_system::prepare()
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
		const auto i = 4 * index + 3;

		const auto a1 = a[i - 3];
		const auto a2 = a[i - 2];
		const auto a3 = a[i - 1];
		const auto a4 = a[i];
		const auto a5 = a[i + 1];
		const auto a6 = a[i + 2];
		const auto a7 = a[i + 3];

		const auto b1 = b[i - 3];
		const auto b2 = b[i - 2];
		const auto b3 = b[i - 1];
		const auto b4 = b[i];
		const auto b5 = b[i + 1];
		const auto b6 = b[i + 2];
		const auto b7 = b[i + 3];

		const auto c1 = c[i - 3];
		const auto c2 = c[i - 2];
		const auto c3 = c[i - 1];
		const auto c4 = c[i];
		const auto c5 = c[i + 1];
		const auto c6 = c[i + 2];
		const auto c7 = c[i + 3];

		const auto c1b0 = c2 * b1;
		const auto c3b2 = c4 * b3;
		const auto c5b4 = c6 * b5;

		const auto b4a3a2 = b5 * a4 * a3;
		const auto c4c3b2 = c5 * c4 * b3;

		const auto b2x = b3 * (a2 * c1 - b2 * b1);
		const auto mb4x = -b5 * (a4 * c3 - b4 * b3);
		const auto b6x = b7 * (a6 * c5 - b6 * b5);

		const auto b2xx = b2x + a3 * c1b0;
		const auto mb4xx = mb4x - a5 * c3b2;
		const auto b6xx = b6x + a7 * c5b4;

		const auto b2xxx = b2xx * c4c3b2;

		ld[index] =
			b6xx * b4a3a2 * b3 * a2 * a1;

		md[index] =
			b6xx * (b2xx * mb4xx + b4a3a2 * c3 * c2 * b1) + b2xxx * b7 * a6 * a5;

		ud[index] = b2xxx * c7 * c6 * b5;

		const auto rho1 = r[i - 3];
		const auto rho2 = r[i - 2];
		const auto rho3 = r[i - 1];
		const auto rho4 = r[i];
		const auto rho5 = r[i + 1];
		const auto rho6 = r[i + 2];
		const auto rho7 = r[i + 3];

		rhs[index + 1] = b6xx * (
			b2xx * (-b5 * (a4 * rho3 - rho4 * b3) - rho5 * c3b2) +
			(b3 * (a2 * rho1 - rho2 * b1) + rho3 * c1b0) * b4a3a2
			) + b2xxx * (b7 * (a6 * rho5 - rho6 * b5) + rho7 * c5b4);
	}

	rhs[0] = df;
	rhs[1] -= ld[0] * df;
	rhs[ec] -= ud.back() * dl;
	rhs[ec + 1] = dl;

	timer.stop();
}

std::vector<number> second_reduced_tridiagonal_system::compute_rest()
{
	auto& full = full_system();
	auto& timer = get_timer();

	timer.start_all_time();

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
		const auto i = 4 * index + 3;

		const auto a1 = a[i - 3];
		const auto a2 = a[i - 2];
		const auto a3 = a[i - 1];

		const auto b1 = b[i - 3];
		const auto b2 = b[i - 2];
		const auto b3 = b[i - 1];

		const auto c1 = c[i - 3];
		const auto c2 = c[i - 2];
		const auto c3 = c[i - 1];

		const auto rho1 = r[i - 3];
		const auto rho2 = r[i - 2];
		const auto rho3 = r[i - 1];

		const auto d0 = rhs[index];;
		const auto d4 = rhs[index + 1];

		const auto c2b1 = b1 * c2;

		const auto d2 = (b3 * (a2 * (rho1 - a1 * d0) - b1 * rho2) +
			c2b1 * (rho3 - c3 * d4)) /
			(b3 * (a2 * c1 - b1 * b2) + c2b1 * a3);

		result[i - 3] = d0;
		result[i - 2] = (rho1 - a1 * d0 - c1 * d2) * (1 / b1);
		result[i - 1] = d2;
		result[i] = (rho3 - a3 * d2 - c3 * d4) * (1 / b3);
	}
	result.back() = dl;

	timer.stop();

	return result;
}
