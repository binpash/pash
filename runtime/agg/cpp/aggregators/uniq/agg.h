#ifndef AGG_UNIQ_H
#define AGG_UNIQ_H

#include "main.h"
#include "util.h"
#include "command_opts.h"
#include <map>
#include <iostream>

constexpr cmd_opts g_options{
    cmd_opt{"count",'c', cmd_opt::Argument::None}
};

static inline void uniq_c() noexcept
{
    size_t padding;
    input1() >> padding;
    padding = input1().tellg();
    input1().seekg(std::ios_base::beg);

    std::string commonKeyCand;
    size_t count2;
    input2() >> count2;
    std::getline(input2(), commonKeyCand);

    input1().seekg(-2, std::ios_base::end);
    while (input1().get() != '\n')
        input1().seekg(-2, std::ios_base::cur);
    
    size_t end_pos = input1().tellg();
    std::string last_line;
    size_t count1;
    input1() >> count1;
    std::getline(input1(), last_line);
    input1().seekg(std::ios_base::beg);

    if (last_line == commonKeyCand)
    {
        stream_copy_n(input1(), output(), end_pos);
        output().width(padding);
        output() << count1 + count2 << commonKeyCand << '\n';
    }
    else
    {
        output() << input1().rdbuf();
        output().width(padding);
        output() << count2 << commonKeyCand << '\n';
    }

    output() << input2().rdbuf();
}
static inline void uniq() noexcept
{    
    std::string commonKeyCand;
    std::getline(input2(), commonKeyCand);
    output() << input1().rdbuf();
    input1().clear();
    input1().seekg(-commonKeyCand.size() - 1, std::ios_base::end);
    std::string last_line;
    std::getline(input1(), last_line);
    if (last_line != commonKeyCand)
        output() << commonKeyCand << '\n';
    output() << input2().rdbuf();
}

void aggregate() noexcept
{
    if (g_options.is_present(0))
    {
        uniq_c();
    }
    else
    {
        uniq();
    }
}

#endif