#include "main.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>

std::ifstream g_in1;
std::ifstream g_in2;
std::ofstream g_out;

int main(int argc, const char* argv[])
{
	if (argc < 4)
	{
		std::cerr << "Usage: " << argv[0] << " input1 input2 output\n";
		return EXIT_FAILURE;
	}

	std::ios_base::sync_with_stdio(false);

	g_in1.open(argv[1], std::ios_base::in  | std::ios_base::binary);
	g_in2.open(argv[2], std::ios_base::in  | std::ios_base::binary);
	g_out.open(argv[3], std::ios_base::out | std::ios_base::binary);

	aggregate();

	return EXIT_SUCCESS;
}

[[nodiscard]] std::istream& input1() noexcept
{
	return g_in1;
}
[[nodiscard]] std::istream& input2() noexcept
{
	return g_in2;
}
[[nodiscard]] std::ostream& output() noexcept
{
	return g_out;
}