#include "main.h"
#include <string>
#include <algorithm>
#include <cstdlib>

inline constexpr cmd_opts g_options{
    cmd_opt{"lines",'l', cmd_opt::Argument::None},
    cmd_opt{"words",'w', cmd_opt::Argument::None},
    cmd_opt{"bytes",'c', cmd_opt::Argument::None},
    cmd_opt{"chars",'m', cmd_opt::Argument::None},
    cmd_opt{"chars",'L', cmd_opt::Argument::None}
};

void aggregate() noexcept
{
    size_t padding;
    input1() >> padding; // dummy read
    padding = input1().tellg();
    input1().seekg(std::ios_base::beg);

    int numbers_to_input = 0;
    for (int i = 0; i < 5; ++i)
        numbers_to_input += g_options.is_present(0); // +1 for each option
    if (numbers_to_input == 0)
        numbers_to_input = 3; // by default there are 3

    for (int i = 0; i < numbers_to_input; ++i)
    {
        size_t count1, count2;
        input1() >> count1;
        input2() >> count2;

        output().width(padding);
        output() << count1 + count2;
        if (i != numbers_to_input - 1)
            output() << ' ';
    }

    output() << '\n';
}
