/*
benchmark.hpp version 20150213
Tancrede Lepoint
Public domain.
*/

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#define BENCHMARK

#ifdef BENCHMARK
#include <chrono>

static auto time_start = std::chrono::system_clock::now();
static auto time_end = std::chrono::system_clock::now();

template <class T>
double get_time_ms(T const &start, T const &end)
{
	auto diff = end - start;
	return (long double)(std::chrono::duration_cast<std::chrono::milliseconds>(diff).count());
}
#define TIME(y, x) time_start = std::chrono::system_clock::now(); \
	x; \
	time_end = std::chrono::system_clock::now(); \
	printf("%s: ", y); \
	std::cout << get_time_ms(time_start, time_end) << "ms" << std::endl;
#else
#define TIME(y, x) x; \
	printf("%s: OK\n", y);
#endif

#endif