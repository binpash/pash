#include "main.h"
#include <string>
#include <algorithm>
#include <cstdlib>

inline constexpr cmd_opts g_options{
    cmd_opt{"lines",'l', cmd_opt::Argument::None},
    cmd_opt{"words",'w', cmd_opt::Argument::None},
    cmd_opt{"bytes",'c', cmd_opt::Argument::None},
    cmd_opt{"chars",'m', cmd_opt::Argument::None},
    cmd_opt{"max-line-length",'L', cmd_opt::Argument::None}
};

void aggregate() noexcept
{
    size_t padding = std::numeric_limits<size_t>::max();
    {
        size_t prev_pos = 0;
        while(true)
        {
            size_t dummy;
            input1() >> dummy;
            if (!input1())
                break;
            size_t new_pos = input1().tellg();
            padding = std::min(new_pos - prev_pos, padding);
            prev_pos = new_pos + 1;
        }
        input1().clear();
        input1().seekg(std::ios_base::beg);
    }

    int numbers_to_input = 0;
    for (int i = 0; i < 4; ++i)
        numbers_to_input += g_options.is_present(i); // +1 for each option
    if (numbers_to_input == 0 && !g_options.is_present(4))
        numbers_to_input = 3; // by default there are 3

    for (int i = 0; i < numbers_to_input; ++i)
    {
        size_t count1, count2;
        input1() >> count1;
        input2() >> count2;

        if (i != 0)
            output() << ' ';
        output().width(padding);
        output() << count1 + count2;
    }

    if (g_options.is_present(4))
    {
        size_t max_len1, max_len2;
        input1() >> max_len1;
        input2() >> max_len2;

        if (numbers_to_input != 0)
            output() << ' ';
        output().width(padding);
        output() << std::max(max_len1, max_len2);
    }

    output() << '\n';
}
