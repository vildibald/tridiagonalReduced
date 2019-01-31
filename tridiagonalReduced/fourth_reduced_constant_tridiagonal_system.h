#pragma once
#include <vector>
#include "number.h"
#include "standard_tridiagonal_system.h"
#include "reduced_tridiagonal_system.h"

class fourth_reduced_constant_tridiagonal_system : public reduced_tridiagonal_system
{
protected:
	reduced_tridiagonal_system* clone_impl() const override;

public:
	explicit fourth_reduced_constant_tridiagonal_system(standard_tridiagonal_system tridiagonal_system);

	
	void prepare() override;

	std::vector<number> compute_rest() override;
};


