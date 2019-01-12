#include "pch.h"
#include "third_reduced_tridiagonal_system.h"
#include "utils.h"

reduced_tridiagonal_system* third_reduced_tridiagonal_system::clone_impl() const
{
	return new third_reduced_tridiagonal_system(*this);
}

third_reduced_tridiagonal_system::third_reduced_tridiagonal_system(
	const standard_tridiagonal_system tridiagonal_system) : reduced_tridiagonal_system(
	std::move(tridiagonal_system), 8)
{
}

void third_reduced_tridiagonal_system::prepare()
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
		const auto i = 8 * index + 7;

		const auto a1 = a[i - 7];
		const auto a2 = a[i - 6];
		const auto a3 = a[i - 5];
		const auto a4 = a[i - 4];
		const auto a5 = a[i - 3];
		const auto a6 = a[i - 2];
		const auto a7 = a[i - 1];
		const auto a8 = a[i];
		const auto a9 = a[i + 1];
		const auto a10 = a[i + 2];
		const auto a11 = a[i + 3];
		const auto a12 = a[i + 4];
		const auto a13 = a[i + 5];
		const auto a14 = a[i + 6];
		const auto a15 = a[i + 7];

		const auto b1 = b[i - 7];
		const auto b2 = b[i - 6];
		const auto b3 = b[i - 5];
		const auto b4 = b[i - 4];
		const auto b5 = b[i - 3];
		const auto b6 = b[i - 2];
		const auto b7 = b[i - 1];
		const auto b8 = b[i];
		const auto b9 = b[i + 1];
		const auto b10 = b[i + 2];
		const auto b11 = b[i + 3];
		const auto b12 = b[i + 4];
		const auto b13 = b[i + 5];
		const auto b14 = b[i + 6];
		const auto b15 = b[i + 7];

		const auto c1 = c[i - 7];
		const auto c2 = c[i - 6];
		const auto c3 = c[i - 5];
		const auto c4 = c[i - 4];
		const auto c5 = c[i - 3];
		const auto c6 = c[i - 2];
		const auto c7 = c[i - 1];
		const auto c8 = c[i];
		const auto c9 = c[i + 1];
		const auto c10 = c[i + 2];
		const auto c11 = c[i + 3];
		const auto c12 = c[i + 4];
		const auto c13 = c[i + 5];
		const auto c14 = c[i + 6];
		const auto c15 = c[i + 7];

		const auto c2b1 = c2 * b1;
		const auto c4b3 = c4 * b3;
		const auto c6b5 = c6 * b5;
		const auto c8b7 = c8 * b7;
		const auto c10b9 = c10 * b9;
		const auto c12b11 = c12 * b11;
		const auto c14b13 = c14 * b13;

		const auto b5a4 = b5 * a4;
		const auto b9a8 = b9 * a8;
		const auto b11a10 = b11 * a10;
		const auto b13a12 = b13 * a12;

		const auto a3c2b1 = a3 * c2b1;
		const auto a5c4b3 = a5 * c4b3;
		const auto a7c6b5 = a7 * c6b5;
		const auto a9c8b7 = a9 * c8b7;
		const auto a11c10b9 = a11 * c10b9;
		const auto a13c12b11 = a13 * c12b11;
		const auto a15c14b13 = a15 * c14b13;

		const auto b9a87 = b9a8 * a7;
		const auto b97a8765 = b9a87 * b7 * a6 * a5;

		const auto b3x = b3 * (a2 * c1 - b2 * b1) + a3c2b1;
		const auto b5x = b5 * (a4 * c3 - b4 * b3) + a5c4b3;
		const auto b7x = b7 * (a6 * c5 - b6 * b5) + a7c6b5;
		const auto b9x = b9 * (a8 * c7 - b8 * b7) + a9c8b7;
		const auto b11x = b11 * (a10 * c9 - b10 * b9) + a11c10b9;
		const auto b13x = b13 * (a12 * c11 - b12 * b11) + a13c12b11;
		const auto b15x = b15 * (a14 * c13 - b14 * b13) + a15c14b13;

		const auto b3y = b3x * c4b3 * c5;
		const auto b7y = b7x * c8b7 * c9;
		const auto b11y = b11x * c12b11 * c13;

		const auto b7xx = -(b7x * (b5a4 * a3c2b1 * c3 - b5x * b3x) + b3y * b7 * a6 * a5);
		const auto b11xx = b11x * (b9a8 * a7c6b5 * c7 - b9x * b7x) + b7y * b11a10 * a9;
		const auto b15xx = -(b15x * (b13a12 * a11c10b9 * c11 - b13x * b11x) + b11y * b15 * a14 * a13
			);

		ld[index] = b15xx * b11x * b7x * b97a8765 * b5a4 * b3 * a3 * a2 * a1;

		md[index] = b15xx * (b11xx * b7xx + b11x * b3y * a7c6b5 * b9a8 * c7 * b7 * a6 * a5) +
			b15x * b7xx * b7x * a11c10b9 * a11c10b9 * a9c8b7 * b13a12 * b11a10 * c11 * c9;

		ud[index] = b11y * b7xx * b7y * c14b13 * c10b9 * c15 * c13;

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

		const auto b7rx = b7 * (a6 * rho5 - rho6 * b5) + c6b5 * rho7;
		const auto b11rx = b11 * (a10 * rho9 - rho10 * b9) + c10b9 * rho11;

		rhs[index + 1] =
			b15xx * (b7xx * (b11x * (b7rx * b9a87 - b7x * (b9 * (a8 * rho7 - rho8 * b7)) +
				b11rx * b7y) + b11x * b97a8765 * (b7x * (b5a4 * a3 * (b3 * (a2 * rho1 - rho2 *
					b1) + c2b1 * rho3) - b3x * (b5 * (a4 * rho3 - rho4 * b3) + c4b3 * rho5)) +
					b3y * b7rx)) + b7xx * b7y * c10b9 * c11 * (b15x * (b11rx * b13a12 * a11 - b11x * (
						b13 * (a12 * rho11 - b11 * rho12) + c12b11 * rho13)) + b11y * (b15 * (a14 * rho13 -
							rho14 * b13) + c14b13 * rho15)));
	}

	rhs[0] = df;
	rhs[1] -= ld[0] * df;
	rhs[ec] -= ud.back() * dl;
	rhs[ec + 1] = dl;

	timer.stop();
}

std::vector<number> third_reduced_tridiagonal_system::compute_rest()
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
		const auto i = 8 * index + 7;

		const auto a1 = a[i - 7];
		const auto a2 = a[i - 6];
		const auto a3 = a[i - 5];
		const auto a4 = a[i - 4];
		const auto a5 = a[i - 3];
		const auto a6 = a[i - 2];
		const auto a7 = a[i - 1];

		const auto b1 = b[i - 7];
		const auto b2 = b[i - 6];
		const auto b3 = b[i - 5];
		const auto b4 = b[i - 4];
		const auto b5 = b[i - 3];
		const auto b6 = b[i - 2];
		const auto b7 = b[i - 1];

		const auto c1 = c[i - 7];
		const auto c2 = c[i - 6];
		const auto c3 = c[i - 5];
		const auto c4 = c[i - 4];
		const auto c5 = c[i - 3];
		const auto c6 = c[i - 2];
		const auto c7 = c[i - 1];

		const auto rho1 = r[i - 7];
		const auto rho2 = r[i - 6];
		const auto rho3 = r[i - 5];
		const auto rho4 = r[i - 4];
		const auto rho5 = r[i - 3];
		const auto rho6 = r[i - 2];
		const auto rho7 = r[i - 1];

		const auto d0 = rhs[index];
		const auto d8 = rhs[index + 1];

		const auto c1b0 = c2 * b1;
		const auto c2b1 = b1 * c2;
		const auto c3b2 = c4 * b3;
		const auto c5b4 = c6 * b5;
		const auto c6b5 = b5 * c6;

		const auto b4a3a2 = b5 * a4 * a3;
		const auto c4c3b2 = c5 * c4 * b3;

		const auto b2x = b3 * (a2 * c1 - b2 * b1);
		const auto mb4x = -b5 * (a4 * c3 - b4 * b3);
		const auto b6x = b7 * (a6 * c5 - b6 * b5);

		const auto b2xx = b2x + a3 * c1b0;
		const auto mb4xx = mb4x - a5 * c3b2;
		const auto b6xx = b6x + a7 * c5b4;

		const auto b2xxx = b2xx * c4c3b2;

		const auto d4 =
			(b6xx * (
				b2xx * (-b5 * (a4 * rho3 - rho4 * b3) - rho5 * c3b2) +
				(b3 * (a2 * rho1 - rho2 * b1) + rho3 * c1b0) * b4a3a2
				) + b2xxx * (b7 * (a6 * rho5 - rho6 * b5) + rho7 * c5b4) -
				b6xx * b4a3a2 * b3 * a2 * a1 * d0 - b2xxx * c7 * c6 * b5 * d8) /
				(b6xx * (b2xx * mb4xx + b4a3a2 * c3 * c2 * b1) + b2xxx * b7 * a6 * a5);

		const auto d2 = (b3 * (a2 * (rho1 - a1 * d0) - b1 * rho2) +
			c2b1 * (rho3 - c3 * d4)) /
			(b3 * (a2 * c1 - b1 * b2) + c2b1 * a3);

		const auto d6 = (b7 * (a6 * (rho5 - a5 * d4) - b5 * rho6) +
			c6b5 * (rho7 - c7 * d8)) /
			(b7 * (a6 * c5 - b5 * b6) + c6b5 * a7);

		result[i - 7] = d0;
		result[i - 6] = (rho1 - a1 * d0 - c1 * d2) * (1 / b1);
		result[i - 5] = d2;
		result[i - 4] = (rho3 - a3 * d2 - c3 * d4) * (1 / b3);
		result[i - 3] = d4;
		result[i - 2] = (rho5 - a5 * d4 - c5 * d6) * (1 / b5);
		result[i - 1] = d6;
		result[i] = (rho7 - a7 * d6 - c7 * d8) * (1 / b7);
	}
	result.back() = dl;

	timer.stop();

	return result;
}