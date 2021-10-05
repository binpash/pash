#ifndef AGG_SORT_H
#define AGG_SORT_H

#include "main.h"
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstring>

inline constexpr cmd_opts g_options{
    // TODO
};

bool compare(const std::string& str1, const std::string& str2)
{
    return str1 < str2;
}

void aggregate() noexcept
{
    std::string input1_top, input2_top;
    bool refill1 = true, refill2 = true;
    while (true)
    {
        if (refill1)
        {
            std::getline(input1(), input1_top);
            if (!input1())
                break;
            refill1 = false;
        }
        if (refill2)
        {
            std::getline(input2(), input2_top);
            if (!input2())
                break;
            refill2 = false;
        }

        if (compare(input1_top, input2_top)) // input1_top < input2_top
        {
            output() << input1_top << '\n';
            refill1 = true;
        }
        else
        {
            output() << input2_top << '\n';
            refill2 = true;
        }
    }

    // At this point exactly one of the inputs is empty
    // and the remaining input from the second one will get outputted
    output() << input1().rdbuf() << input2().rdbuf();
}

#endif // AGG_SORT_H