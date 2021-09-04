#ifndef AGG_TR_H
#define AGG_TR_H

#include "main.h"
#include "util.h"
#include "command_opts.h"

constexpr cmd_opts g_options{
    cmd_opt{"complement",'c', cmd_opt::Argument::None},
    cmd_opt{"truncate-set1",'t', cmd_opt::Argument::None}, // linux exclusive
    cmd_opt{"delete",'d', cmd_opt::Argument::None},
    cmd_opt{"squeeze-repeats",'s', cmd_opt::Argument::None},
    cmd_opt{"",'u', cmd_opt::Argument::None} // bsd exclusive
};

inline static void tr() noexcept
{
    nyi_error("The -s option isn't supported yet.");
}

void aggregate() noexcept
{
    if (!g_options.is_present(3))
        output() << input1().rdbuf() << input2().rdbuf();
    else
        tr();
}

#endif
