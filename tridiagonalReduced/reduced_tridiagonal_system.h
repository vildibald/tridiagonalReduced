#pragma once
#include "tridiagonal_system.h"
#include "standard_tridiagonal_system.h"

class reduced_tridiagonal_system : public tridiagonal_system
{
	standard_tridiagonal_system full_system_;

protected:

	reduced_tridiagonal_system(standard_tridiagonal_system full_system, int factor);

	virtual reduced_tridiagonal_system* clone_impl() const = 0;

public:

	standard_tridiagonal_system& full_system();

	virtual ~reduced_tridiagonal_system();

	auto clone() const { return std::unique_ptr<reduced_tridiagonal_system>(clone_impl()); };

	int reduction_factor() const;

	virtual void prepare() = 0;

	virtual std::vector<number> compute_rest() = 0;

	std::vector<number> solve() override;
};
