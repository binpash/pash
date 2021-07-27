#include "main.h"
#include "agg.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <tuple>

[[noreturn]] void nyi_error(const char* message) noexcept
{
	std::cerr << "Not yet implemented error: " << message << '\n';
	exit(EXIT_FAILURE);
}

std::pair<int, char**> parse_options(int argc, char** argv, const opt_holder& options) noexcept
{	
	int ch, longind = 0;
	while ((ch = getopt_long(argc, argv, options.optstring, options.long_options, &longind)) != -1)
	{
		size_t opt_idx;
		switch(ch)
		{
		case 0:
			opt_idx = (options.long_options[longind].flag - options.present);	
			break;
		case '?':
			nyi_error(("Unsupported flag " + std::to_string(ch)).c_str());
		default:
			opt_idx = options.map[ch];
			break;
		}
		if (optarg)
			options.args[opt_idx] = optarg;
		options.present[opt_idx] = 1;
	}
	return {argc - optind, argv + optind};
}

// Written to by main
// Read by functions below it
std::ifstream g_in1;
std::ifstream g_in2;
std::ofstream g_out;
int g_argc;
char** g_argv;

int main(int _argc, char** _argv)
{
	std::ios_base::sync_with_stdio(false);

	if (_argc < 4)
	{
		std::cerr << "Usage: " << _argv[0] << " input1 input2 output [OPTION]...\n";
		return EXIT_FAILURE;
	}

	g_in1.open(_argv[1], std::ios_base::in  | std::ios_base::binary);
	g_in2.open(_argv[2], std::ios_base::in  | std::ios_base::binary);
	g_out.open(_argv[3], std::ios_base::out | std::ios_base::binary);

	if (g_in1.peek() == EOF)
	{
		output() << input2().rdbuf();
		return EXIT_SUCCESS;
	}
	if (g_in2.peek() == EOF)
	{
		output() << input1().rdbuf();
		return EXIT_SUCCESS;
	}

	std::tie(g_argc, g_argv) = parse_options(_argc - 3, _argv + 3, g_options);

	aggregate();

	return EXIT_SUCCESS;
}

[[nodiscard]] size_t argc() noexcept
{
	return g_argc;
}
[[nodiscard]] char** argv() noexcept
{
	return g_argv;
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
