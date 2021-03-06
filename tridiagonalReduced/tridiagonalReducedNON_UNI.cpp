// tridiagonalReduced.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "number.h"
#include <vector>
#include "comparison_benchmark_result.h"
#include <functional>
#include "vector_extensions.h"
#include "standard_tridiagonal_system.h"
#include "first_reduced_tridiagonal_system.h"
#include <algorithm>
#include "second_reduced_tridiagonal_system.h"
#include "third_reduced_tridiagonal_system.h"
#include <random>
#include "first_reduced_constant_tridiagonal_system.h"
#include "second_reduced_constant_tridiagonal_system.h"
#include "third_reduced_constant_tridiagonal_system.h"
#include "fourth_reduced_constant_tridiagonal_system.h"

int round_num_knots(const int num_knots)
{
	const auto wanted = 64;
	const auto remainder = num_knots % wanted;
	if (remainder != wanted - 1)
	{
		return (1 + num_knots / wanted) * wanted - 1;
	}
	return num_knots;
}

std::vector<number> test_vector(const number from, const number to, const int size)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);

	std::vector<number> vector(size);
	for (int i = 0; i < size; ++i)
	{
		vector[i] = dist(mt);
	}
	return vector;
}

comparison_benchmark_result benchmark_general_solver(const int num_iterations, int num_knots)
{
	num_knots = round_num_knots(num_knots);

	auto main_diagonal = test_vector(-1000, 1000, num_knots);
	auto lower_diagonal = ve::map(main_diagonal, ve::quarterifier());
	auto upper_diagonal = ve::map(main_diagonal, ve::eightifier());
	auto right_side = test_vector(-100, 100, num_knots);
	const auto d_first = -0.14;
	const auto d_last = 0.14;

	auto full = standard_tridiagonal_system(std::move(lower_diagonal),
	                                        std::move(main_diagonal), std::move(upper_diagonal),
	                                        std::move(right_side), d_first, d_last);
	auto first_reduced = first_reduced_tridiagonal_system(full);
	auto second_reduced = second_reduced_tridiagonal_system(full);
	auto third_reduced = third_reduced_tridiagonal_system(full);

	std::vector<number> calculated_results;
	long long full_time = 0;
	long long first_reduced_time = 0;
	long long second_reduced_time = 0;
	long long third_reduced_time = 0;

	calculated_results.reserve(num_iterations * 6);

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = full;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		full_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = first_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		first_reduced_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = second_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		second_reduced_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = third_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		third_reduced_time += solver.get_timer().execution_time();
	}

	std::cout << "Ignore " << calculated_results[0] << std::endl;
	comparison_benchmark_result benchmark_result;
	benchmark_result
		.add(full_time)
		.add(first_reduced_time)
		.add(second_reduced_time)
		.add(third_reduced_time);
	return benchmark_result;
}

comparison_benchmark_result benchmark_constant_solver(const int num_iterations, int num_knots)
{
	num_knots = round_num_knots(num_knots);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 400.0);

	auto main_diagonal = std::vector<number>(num_knots, dist(mt));
	auto lower_diagonal = ve::map(main_diagonal, ve::quarterifier());
	auto upper_diagonal = ve::map(main_diagonal, ve::eightifier());
	auto right_side = ve::map(main_diagonal, ve::twicier());
	const auto d_first = -0.14;
	const auto d_last = 0.14;

	auto full = standard_tridiagonal_system(std::move(lower_diagonal),
		std::move(main_diagonal), std::move(upper_diagonal),
		std::move(right_side), d_first, d_last);
	auto first_reduced = first_reduced_constant_tridiagonal_system(full);
	auto second_reduced = second_reduced_constant_tridiagonal_system(full);
	auto third_reduced = third_reduced_constant_tridiagonal_system(full);
	auto fourth_reduced = fourth_reduced_constant_tridiagonal_system(full);

	std::vector<number> calculated_results;
	long long full_time = 0;
	long long first_reduced_time = 0;
	long long second_reduced_time = 0;
	long long third_reduced_time = 0;
	long long fourth_reduced_time = 0;

	calculated_results.reserve(num_iterations * 6);

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = full;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		full_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = first_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		first_reduced_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = second_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		second_reduced_time += solver.get_timer().execution_time();
	}

	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = third_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		third_reduced_time += solver.get_timer().execution_time();
	}


	for (int i = 0; i < num_iterations; i++)
	{
		auto solver = fourth_reduced;
		auto result = solver.solve();
		calculated_results.push_back(result[1]);
		fourth_reduced_time += solver.get_timer().execution_time();
	}

	std::cout << "Ignore " << calculated_results[0] << std::endl;
	comparison_benchmark_result benchmark_result;
	benchmark_result
		.add(full_time)
		.add(first_reduced_time)
		.add(second_reduced_time)
		.add(third_reduced_time);
	return benchmark_result;
}

void print_general_matrix_result(comparison_benchmark_result& result)
{
	std::cout << "0th Reduced : " << result[0] << std::endl;
	std::cout << "1st Reduced : " << result[1] << std::endl;
	std::cout << "2nd Reduced : " << result[2] << std::endl;
	std::cout << "3rd Reduced : " << result[3] << std::endl;
	std::cout << ".........." << std::endl;
	std::cout << "Speedup 1st Reduced: " << result.ratio(0, 1) << std::endl;
	std::cout << "Speedup 2nd Reduced: " << result.ratio(0, 2) << std::endl;
	std::cout << "Speedup 3rd Reduced: " << result.ratio(0, 3) << std::endl;
}

void print_constant_matrix_result(comparison_benchmark_result& result)
{
	std::cout << "0th Reduced : " << result[0] << std::endl;
	std::cout << "1st Reduced : " << result[1] << std::endl;
	std::cout << "2nd Reduced : " << result[2] << std::endl;
	std::cout << "3rd Reduced : " << result[3] << std::endl;
	std::cout << "4th Reduced : " << result[4] << std::endl;
	std::cout << ".........." << std::endl;
	std::cout << "Speedup 1st Reduced: " << result.ratio(0, 1) << std::endl;
	std::cout << "Speedup 2nd Reduced: " << result.ratio(0, 2) << std::endl;
	std::cout << "Speedup 3rd Reduced: " << result.ratio(0, 3) << std::endl;
	std::cout << "Speedup 4th Reduced: " << result.ratio(0, 4) << std::endl;
}

void perform_general_solver_benchmark()
{
	int num_iterations;
	int num_knots;
	std::cout << "Tridiagonal solver benchmark" << std::endl << std::endl;
	std::cout << "Enter number of iterations: " << std::endl;
	std::cin >> num_iterations;
	std::cout << "Enter number of knots: " << std::endl;
	std::cin >> num_knots;
	std::cin.get();

	auto result = benchmark_general_solver(num_iterations, num_knots);
	print_general_matrix_result(result);
}

void perform_constant_solver_benchmark()
{
	int num_iterations;
	int num_knots;
	std::cout << "Tridiagonal solver benchmark" << std::endl << std::endl;
	std::cout << "Enter number of iterations: " << std::endl;
	std::cin >> num_iterations;
	std::cout << "Enter number of knots: " << std::endl;
	std::cin >> num_knots;
	std::cin.get();

	auto result = benchmark_constant_solver(num_iterations, num_knots);
	print_constant_matrix_result(result);
}

void general_matrix_equality_comparison(int num_knots, const number range)
{
	num_knots = round_num_knots(num_knots);

	auto main_diagonal = test_vector(-range, range, num_knots);
	auto lower_diagonal = ve::map(main_diagonal, ve::quarterifier());
	auto upper_diagonal = ve::map(main_diagonal, ve::eightifier());
	auto right_side = std::vector<number>(num_knots, 1);
	const auto d_first = -0.14;
	const auto d_last = 0.14;

	auto full = standard_tridiagonal_system(std::move(lower_diagonal),
	                                        std::move(main_diagonal), std::move(upper_diagonal),
	                                        std::move(right_side), d_first, d_last);
	auto first_reduced = first_reduced_tridiagonal_system(full);
	auto second_reduced = second_reduced_tridiagonal_system(full);
	auto third_reduced = third_reduced_tridiagonal_system(full);

	const auto result_full = full.solve();
	const auto result_first_reduced = first_reduced.solve();
	const auto result_second_reduced = second_reduced.solve();
	auto result_third_reduced = third_reduced.solve();

	auto min_diff_r1 = (std::numeric_limits<number>::max)(),
		max_diff_r1 = (std::numeric_limits<number>::min)(),
		min_diff_r2 = (std::numeric_limits<number>::max)(),
		max_diff_r2 = (std::numeric_limits<number>::min)(),
		min_diff_r3 = (std::numeric_limits<number>::max)(),
		max_diff_r3 = (std::numeric_limits<number>::min)();

	for (int i = 0; i < result_full.size(); ++i)
	{
		if (result_third_reduced[i] == std::numeric_limits<number>::quiet_NaN())
		{
			result_third_reduced[i] = (std::numeric_limits<number>::max)() * 0.5;
		}
		min_diff_r1 = (std::min)(min_diff_r1,
		                         std::abs(result_full[i] - result_first_reduced[i]));
		max_diff_r1 = (std::max)(max_diff_r1,
		                         std::abs(result_full[i] - result_first_reduced[i]));

		min_diff_r2 = (std::min)(min_diff_r2,
		                         std::abs(result_full[i] - result_second_reduced[i]));
		max_diff_r2 = (std::max)(max_diff_r2,
		                         std::abs(result_full[i] - result_second_reduced[i]));

		min_diff_r3 = (std::min)(min_diff_r3,
		                         std::abs(result_full[i] - result_third_reduced[i]));
		max_diff_r3 = (std::max)(max_diff_r3,
		                         std::abs(result_full[i] - result_third_reduced[i]));

		if (num_knots < 1000) {
			std::cout << "\nColumn " << i << ":\n"
				<< "0th R.: \t" << result_full[i] << '\n';
			std::cout << "1st R.: \t" << result_first_reduced[i] << '\n';
			std::cout << "2nd R.: \t" << result_second_reduced[i] << '\n';
			std::cout << "3rd R.: \t" << result_third_reduced[i] << '\n';
		}
	}
	std::cout << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 1st R.: " << min_diff_r1 << std::endl;
	std::cout << "Max diff 1st R.: " << max_diff_r1 << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 2nd R.: " << min_diff_r2 << std::endl;
	std::cout << "Max diff 2nd R.: " << max_diff_r2 << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 3rd R.: " << min_diff_r3 << std::endl;
	std::cout << "Max diff 3rd R.: " << max_diff_r3 << std::endl;
}

void constant_matrix_equality_comparison(int num_knots)
{
	num_knots = round_num_knots(num_knots);
	auto main_diagonal = std::vector<number>(num_knots, 432);
	auto lower_diagonal = std::vector<number>(num_knots, 123);
	auto upper_diagonal = std::vector<number>(num_knots, 187);
	auto right_side = std::vector<number>(num_knots, 1);
	const auto d_first = -0.14;
	const auto d_last = 0.14;

	auto full = standard_tridiagonal_system(std::move(lower_diagonal),
		std::move(main_diagonal), std::move(upper_diagonal),
		std::move(right_side), d_first, d_last);
	auto first_reduced = first_reduced_constant_tridiagonal_system(full);
	auto second_reduced = second_reduced_constant_tridiagonal_system(full);
	auto third_reduced = third_reduced_constant_tridiagonal_system(full);
	auto fourth_reduced = fourth_reduced_constant_tridiagonal_system(full);

	const auto result_full = full.solve();
	const auto result_first_reduced = first_reduced.solve();
	const auto result_second_reduced = second_reduced.solve();
	const auto result_third_reduced = third_reduced.solve();
	auto result_fourth_reduced = fourth_reduced.solve();

	auto min_diff_r1 = (std::numeric_limits<number>::max)(),
		max_diff_r1 = (std::numeric_limits<number>::min)(),
		min_diff_r2 = (std::numeric_limits<number>::max)(),
		max_diff_r2 = (std::numeric_limits<number>::min)(),
		min_diff_r3 = (std::numeric_limits<number>::max)(),
		max_diff_r3 = (std::numeric_limits<number>::min)(),
		min_diff_r4 = (std::numeric_limits<number>::max)(),
		max_diff_r4 = (std::numeric_limits<number>::min)();

	for (int i = 0; i < result_full.size(); ++i)
	{
		if(result_fourth_reduced[i] == std::numeric_limits<number>::quiet_NaN())
		{
			result_fourth_reduced[i] = (std::numeric_limits<number>::max)() * 0.5;
		}
		min_diff_r1 = (std::min)(min_diff_r1,
			std::abs(result_full[i] - result_first_reduced[i]));
		max_diff_r1 = (std::max)(max_diff_r1,
			std::abs(result_full[i] - result_first_reduced[i]));

		min_diff_r2 = (std::min)(min_diff_r2,
			std::abs(result_full[i] - result_second_reduced[i]));
		max_diff_r2 = (std::max)(max_diff_r2,
			std::abs(result_full[i] - result_second_reduced[i]));

		min_diff_r3 = (std::min)(min_diff_r3,
			std::abs(result_full[i] - result_third_reduced[i]));
		max_diff_r3 = (std::max)(max_diff_r3,
			std::abs(result_full[i] - result_third_reduced[i]));

		min_diff_r4 = (std::min)(min_diff_r4,
			std::abs(result_full[i] - result_fourth_reduced[i]));
		max_diff_r4 = (std::max)(max_diff_r4,
			std::abs(result_full[i] - result_fourth_reduced[i]));

		if (num_knots < 1000) {
			std::cout << "\nColumn " << i << ":\n"
				<< "0th R.: \t" << result_full[i] << '\n';
			std::cout << "1st R.: \t" << result_first_reduced[i] << '\n';
			std::cout << "2nd R.: \t" << result_second_reduced[i] << '\n';
			std::cout << "3rd R.: \t" << result_third_reduced[i] << '\n';
			std::cout << "4th R.: \t" << result_fourth_reduced[i] << '\n';
		}
	}
	std::cout << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 1st R.: " << min_diff_r1 << std::endl;
	std::cout << "Max diff 1st R.: " << max_diff_r1 << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 2nd R.: " << min_diff_r2 << std::endl;
	std::cout << "Max diff 2nd R.: " << max_diff_r2 << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 3rd R.: " << min_diff_r3 << std::endl;
	std::cout << "Max diff 3rd R.: " << max_diff_r3 << std::endl;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Min diff 4th R.: " << min_diff_r4 << std::endl;
	std::cout << "Max diff 4th R.: " << max_diff_r4 << std::endl;
}

void perform_general_equality_comparison()
{
	int num_knots;
	number range;
	std::cout << "Equality comparison" << std::endl << std::endl;
	std::cout << "Enter number of knots: " << std::endl;
	std::cin >> num_knots;
	std::cout << "Enter max absolute size of the coefficients (decimal): " << std::endl;
	std::cin >> range;
	std::cin.get();

	general_matrix_equality_comparison(num_knots, range);
}

void perform_constant_equality_comparison()
{
	int num_knots;
	std::cout << "Equality comparison" << std::endl << std::endl;
	std::cout << "Enter number of knots: " << std::endl;
	std::cin >> num_knots;
	std::cin.get();

	constant_matrix_equality_comparison(num_knots);
}

int main()
{
	while (true)
	{
		//std::cout << clock();
		// Console clear ...
		// ... for Windows,
		system("cls");
		// ... for Linux/Unix.
		//system("clear");
		std::cout << "1: General tridiagonal solver benchmark." << std::endl;
		std::cout << "2: General tridiagonal solver equality comparison." << std::endl;
		std::cout << "3: Constant tridiagonal solver benchmark." << std::endl;
		std::cout << "4: Constant tridiagonal solver equality comparison." << std::endl;
		std::cout << "Q: End program" << std::endl;
		char input;
		std::cin >> input;
		std::cin.get();
		std::cout << std::endl << "---------------" << std::endl;

		switch (input)
		{
		case '1':
			perform_general_solver_benchmark();
			break;
		case '2':
			perform_general_equality_comparison();
			break;
		case '3':
			perform_constant_solver_benchmark();
			break;
		case '4':
			perform_constant_equality_comparison();
			break;
		case 'q':
		case 'Q':
			return 0;
		default:;
		}

		std::cout << "===================" << std::endl;

		system("pause");
	}
}