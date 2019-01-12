#pragma once
#include <vector>
#include <functional>

namespace ve
{
	template <typename S, typename R>
	std::vector<R> map(std::vector<S>& source, std::function<R(S)>& operation)
	{
		std::vector<R> result(source.size());
		for (int i = 0; i < source.size(); ++i)
		{
			result[i] = operation(source[i]);
		}
		return result;
	}

	template <typename S, typename R>
	std::vector<R> map(std::vector<S>& source, std::function<R(S)>&& operation)
	{
		std::vector<R> result(source.size());
		for (int i = 0; i < source.size(); ++i)
		{
			result[i] = operation(source[i]);
		}
		return result;
	}

	std::function<number(number)> quarterifier()
	{
		return [](const number element)
		{
			return element / 4.0;
		};
	}

	std::function<number(number)> eightifier()
	{
		return [](const number element) { return element / 8.0; };
	}

	std::function<number(number)> twicier()
	{
		return [](const number element) { return element * 2; };
	}
}
