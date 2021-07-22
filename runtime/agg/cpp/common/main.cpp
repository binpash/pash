#include "main.h"
#include "agg.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

std::ifstream g_in1;
std::ifstream g_in2;
std::ofstream g_out;

void parse_options(int argc, char** argv, const opt_holder& options) noexcept
{
	int ch, longind;
	while ((ch = getopt_long(argc, argv, options.optstring, options.long_options, &longind)) != -1)
	{
		switch(ch)
		{
		case 0:
		{
			auto* opt = reinterpret_cast<cmd_opt*>(reinterpret_cast<char*>(options.long_options[longind].flag) - offsetof(cmd_opt, present));
			size_t optarglen = strlen(optarg);
			opt->arg = new char[optarglen + 1];
			strcpy(opt->arg, optarg);
			break;
		}
		case '?':
			nyi_error(("Unsupported flag " + std::to_string(ch)).c_str());
		default:
			break;
		}
	}
}

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		std::cerr << "Usage: " << argv[0] << " input1 input2 output [OPTION]...\n";
		return EXIT_FAILURE;
	}

	parse_options(argc - 4, argv + 4, g_options);

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

[[noreturn]] void nyi_error(const char* message) noexcept
{
	std::cerr << "Not yet implemented error: " << message << '\n';
	exit(EXIT_FAILURE);
}