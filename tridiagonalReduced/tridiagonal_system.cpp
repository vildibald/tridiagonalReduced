#include "pch.h"
#include "tridiagonal_system.h"

tridiagonal_system::~tridiagonal_system() = default;

tridiagonal_system::tridiagonal_system(const std::vector<number>& lower_diagonal,
                                       const std::vector<number>& main_diagonal,
                                       const std::vector<number>& upper_diagonal,
                                       const std::vector<number>& right_side, number d_first, number d_last,
                                       const timer& timer): lower_diagonal_(lower_diagonal),
                                                            main_diagonal_(main_diagonal),
                                                            upper_diagonal_(upper_diagonal),
                                                            right_side_(right_side),
                                                            d_first_(d_first),
                                                            d_last_(d_last),
                                                            timer_(timer)
{
}

tridiagonal_system::tridiagonal_system(
	const std::vector<number> lower_diagonal,
	const std::vector<number> main_diagonal,
	const std::vector<number> upper_diagonal,
	const std::vector<number> right_side,
	const number d_first, const number d_last) : tridiagonal_system(
		std::move(lower_diagonal),
		std::move(main_diagonal),
		std::move(upper_diagonal),
		std::move(right_side),
		d_first,
		d_last,
		timer())
{
}

int tridiagonal_system::equation_count() const
{
	return main_diagonal_.size();
}

std::vector<number>& tridiagonal_system::lower_diagonal()
{
	return lower_diagonal_;
}

std::vector<number>& tridiagonal_system::main_diagonal()
{
	return main_diagonal_;
}

std::vector<number>& tridiagonal_system::upper_diagonal()
{
	return upper_diagonal_;
}

std::vector<number>& tridiagonal_system::right_side()
{
	return right_side_;
}

number tridiagonal_system::d_first() const
{
	return d_first_;
}

number tridiagonal_system::d_last() const
{
	return d_last_;
}

timer& tridiagonal_system::get_timer()
{
	return timer_;
}