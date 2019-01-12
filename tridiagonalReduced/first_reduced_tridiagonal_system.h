#pragma once
#include <vector>
#include "number.h"
#include "standard_tridiagonal_system.h"
#include "reduced_tridiagonal_system.h"

class first_reduced_tridiagonal_system final : public reduced_tridiagonal_system
{
protected:
	reduced_tridiagonal_system* clone_impl() const override;

public:
	explicit first_reduced_tridiagonal_system(standard_tridiagonal_system tridiagonal_system);

	void prepare() override;

	std::vector<number> compute_rest() override;
};
