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
    bool invalid1 = true, invalid2 = true;
    while (true)
    {
        if (invalid1)
        {
            std::getline(input1(), input1_top);
            if (!input1())
                break;
            invalid1 = false;
        }
        if (invalid2)
        {
            std::getline(input2(), input2_top);
            if (!input2())
                break;
            invalid2 = false;
        }

        if (compare(input1_top, input2_top)) // input1_top < input2_top
        {
            output() << input1_top << '\n';
            invalid1 = true;
        }
        else
        {
            output() << input2_top << '\n';
            invalid2 = true;
        }
    }

    if(!invalid1)
        output() << input1_top << '\n';
    if(!invalid2)
        output() << input2_top << '\n';

    // At this point at least one of the inputs
    // is empty, so the other one can sefely be forwarded
    output() << input1().rdbuf() << input2().rdbuf();
}

#endif // AGG_SORT_H