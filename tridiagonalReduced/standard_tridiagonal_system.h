#pragma once
#include <vector>
#include "number.h"
#include "timer.h"
#include "tridiagonal_system.h"

class standard_tridiagonal_system final : public tridiagonal_system
{
public:
	standard_tridiagonal_system(
		std::vector<number> lower_diagonal,
		std::vector<number> main_diagonal,
		std::vector<number> upper_diagonal,
		std::vector<number> right_side,
		number d_first,
		number d_last,
		timer timer);

	standard_tridiagonal_system(
		std::vector<number> lower_diagonal,
		std::vector<number> main_diagonal,
		std::vector<number> upper_diagonal,
		std::vector<number> right_side,
		number d_first,
		number d_last);

	std::vector<number> solve() override;
};
