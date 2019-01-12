#pragma once
#include "number.h"
#include <vector>
#include "timer.h"

class equation_system
{
public:
	virtual ~equation_system() = default;

	virtual std::vector<number> solve() = 0;

	virtual timer& get_timer() = 0;
};

