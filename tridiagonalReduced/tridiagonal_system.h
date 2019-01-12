#pragma once
#include "timer.h"
#include "number.h"
#include <vector>
#include "equation_system.h"
#include <memory>

class tridiagonal_system : public equation_system
{
	std::vector<number> lower_diagonal_;
	std::vector<number> main_diagonal_;
	std::vector<number> upper_diagonal_;
	std::vector<number> right_side_;
	number d_first_;
	number d_last_;
	timer timer_;

protected:
	tridiagonal_system(const std::vector<number>& lower_diagonal,
	                   const std::vector<number>& main_diagonal,
	                   const std::vector<number>& upper_diagonal,
	                   const std::vector<number>& right_side, number d_first, number d_last,
	                   const timer& timer);

	tridiagonal_system(std::vector<number> lower_diagonal,
	                   std::vector<number> main_diagonal,
	                   std::vector<number> upper_diagonal,
	                   std::vector<number> right_side,
	                   number d_first,
	                   number d_last);

	

public:
	virtual ~tridiagonal_system();

	int equation_count() const;

	std::vector<number>& lower_diagonal();

	std::vector<number>& main_diagonal();

	std::vector<number>& upper_diagonal();

	std::vector<number>& right_side();

	number d_first() const;

	number d_last() const;

	timer& get_timer() override;
};
