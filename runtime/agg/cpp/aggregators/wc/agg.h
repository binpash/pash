#include "main.h"
#include <string>
#include <algorithm>
#include <cstdlib>

inline constexpr cmd_opts g_options{};

void aggregate() noexcept
{
    size_t line_count1, word_count1, char_count1;
    input1() >> line_count1;
    size_t padding = input1().tellg();
    input1() >> word_count1 >> char_count1;

    size_t line_count2, word_count2, char_count2;
    input2() >> line_count2 >> word_count2 >> char_count2;

    output().width(padding);
    output() << line_count1 + line_count2 << ' ';
    output().width(padding);
    output() << word_count1 + word_count2 << ' ';
    output().width(padding);
    output() << char_count1 + char_count2 << '\n';
}
