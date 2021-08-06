#ifndef AGG_GREP_H
#define AGG_GREP_H

#include "main.h"
#include "util.h"
#include "command_opts.h"

constexpr cmd_opts g_options{
    // These options require special implementations
    cmd_opt{"count", 'c', cmd_opt::Argument::None},
    //cmd_opt{"max-count", 'm', cmd_opt::Argument::Required},
    //cmd_opt{"byte-offset", 'b', cmd_opt::Argument::None},
    //cmd_opt{"line-number", 'n', cmd_opt::Argument::None},
    //cmd_opt{"label", '\0', cmd_opt::Argument::Required},

    // These options have no effect, they're here only not to
    // throw an unrecognized option error
    cmd_opt{"extended-regexp", 'E', cmd_opt::Argument::None},
    cmd_opt{"fixed-strings", 'F', cmd_opt::Argument::None},
    cmd_opt{"basic-regexp", 'G', cmd_opt::Argument::None},
    cmd_opt{"perl-regexp", 'P', cmd_opt::Argument::None},
    cmd_opt{"regexp", 'e', cmd_opt::Argument::Required},
    cmd_opt{"file", 'f', cmd_opt::Argument::Required},
    cmd_opt{"ignore-case", 'i', cmd_opt::Argument::None},
    cmd_opt{"no-ignore-case", '\0', cmd_opt::Argument::None},
    cmd_opt{"word-regexp", 'w', cmd_opt::Argument::None},
    cmd_opt{"line-regexp", 'x', cmd_opt::Argument::None},
    cmd_opt{"no-messages", 's', cmd_opt::Argument::None},
    cmd_opt{"invert-match", 'v', cmd_opt::Argument::None},
    cmd_opt{"line-buffered", '\0', cmd_opt::Argument::None},
    cmd_opt{"only-matching", 'o', cmd_opt::Argument::None},
    cmd_opt{"quiet", 'q', cmd_opt::Argument::None},
    cmd_opt{"silent", '\0', cmd_opt::Argument::None},
    cmd_opt{"binary-files", '\0', cmd_opt::Argument::Required},
    cmd_opt{"text", 'a', cmd_opt::Argument::None},
    cmd_opt{"",'I', cmd_opt::Argument::None}
};

void aggregate() noexcept
{
    if (!g_options.is_present(0))
        output() << input1().rdbuf() << input2().rdbuf();
    else
    {
        size_t count1, count2;
        input1() >> count1;
        input2() >> count2;
        output() << count1 + count2 << '\n';
    }
}

#endif
