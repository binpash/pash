#include "main.h"
#include <algorithm>
#include <iostream>
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

[[nodiscard]] std::optional<std::string> input(std::ifstream& in) noexcept
{
    if (!in)
        return std::nullopt;
    
    std::string s;
    std::getline(in, s)

    return s;
}
[[nodiscard]] std::string inputAll(std::ifstream& in) noexcept
{
    if (!in)
        return {};
    
    std::string contents;
	in.seekg(0, std::ios::end);
	contents.resize(in.tellg());
	in.seekg(0, std::ios::beg);
	in.read(&contents[0], contents.size());

	return contents;
}
[[nodiscard]] std::pair<std::string, std::vector<std::string_view>> inputAllLines(std::ifstream& in) noexcept
{
    std::string contents{inputAll(in)};
    int lineCount = std::count(contents.begin(), contents.end(), '\n');
    std::vector<std::string_view> lines{};
    lines.reserve(lineCount);
    for (int i = 0, first = 0; i < contents.size(); ++i)
    {
        if (contents[i] == '\n')
        {
            lines.emplace_back(contents + first, i - first);
            first = i + 1;
        }
    }
    return {contents, lines};
}

[[nodiscard]] std::optional<std::string> input1() noexcept { return input(g_in1); }
[[nodiscard]] std::optional<std::string> input2() noexcept { return input(g_in2); }
[[nodiscard]] std::string input1All() noexcept { return inputAll(g_in1); }
[[nodiscard]] std::string input2All() noexcept { return inputAll(g_in2); }

[[nodiscard]] std::pair<std::string, std::vector<std::string_view>> input1AllLines() noexcept { return inputAllLines(g_in1); }
[[nodiscard]] std::pair<std::string, std::vector<std::string_view>> input2AllLines() noexcept { return inputAllLines(g_in2); }
void output(const std::string& s) noexcept { output(s.c_str(), s.size()); }
void output(const char* s, size_t len) noexcept
{
    g_out.write(s, len);
    if (s[len - 1] != '\n')
        g_out.put('\n');
}