#ifndef COMMAND_OPTS_H
#define COMMAND_OPTS_H

#include <vector>
#include <array>
#include <initializer_list>
#include <getopt.h>

struct opt_holder
{
    option* long_options;
    const char* optstring;
};

struct cmd_opt
{
    constexpr cmd_opt(const char* long_name, char abbreviation) noexcept:
        long_name{long_name},
        abbreviation{abbreviation},
        present{false},
        arg{nullptr}
    {
    }

    const char* long_name;
    char abbreviation;
    int present;
    char* arg;
};

template <size_t Size>
class cmd_opts : public opt_holder, public std::array<cmd_opt, Size>
{
public:
    constexpr cmd_opts(std::initializer_list<cmd_opt> ar) noexcept :
        opt_holder{long_opts.data(), optstring.data()},
        std::array<cmd_opt, Size>{il_to_array<Size>(ar)},
        long_opts{[&](){
            std::array<option, Size + 1> long_opts{};
            size_t current_idx = 0;
            auto& base = static_cast<std::array<cmd_opt, Size>&>(*this);
            for (auto& opt : base)
            {
                if (opt.long_name)
                {
                    long_opts[current_idx++] = {
                        opt.long_name,
                        optional_argument,
                        &opt.present,
                        1
                    };
                }
            }
            return long_opts;
        }()},
        optstring{[&](){
            std::array<char, 3 * Size + 1> optstring{};
            size_t current_idx = 0;
            auto& base = static_cast<std::array<cmd_opt, Size>&>(*this);
            for (auto& opt : base)
            {
                if (opt.abbreviation)
                {
                    optstring[current_idx++] = opt.abbreviation;
                    optstring[current_idx++] = ':';
                    optstring[current_idx++] = ':';
                }
            }
            return optstring;
        }()}
        {
        }

        // Technically, this should be uncommented to not have a memory leak
        // however, then the code won't compile as the destructor is non-trivial
        // Fortunately all these allocations will be freed by the OS anyway on exit
        // so this shouldn't be a problem.
        //
        // ~cmd_opts()
        // {
        //     auto& base = static_cast<std::array<cmd_opt, Size>&>(*this);
        //     for (auto& opt : base)
        //     {
        //         delete[] opt.arg; // safe if nullptr
        //         opt.arg = nullptr;
        //     }
        // }
private:
    std::array<option, Size + 1> long_opts;
    std::array<char, 3 * Size + 1> optstring;

    template<size_t N, typename T, std::size_t... I>
    [[nodiscard]] constexpr auto il_to_array_impl(const std::initializer_list<T>& il, std::index_sequence<I...>) noexcept
    {
        return std::array<T, N>{std::data(il)[I]...};
    }
    template<size_t N, typename T>
    [[nodiscard]] constexpr auto il_to_array(const std::initializer_list<T>& il) noexcept
    {
        return il_to_array_impl<N>(il, std::make_index_sequence<N>{});
    }
};

template <typename... T>
cmd_opts(T...) -> cmd_opts<sizeof...(T)>;
cmd_opts() -> cmd_opts<0>;

#endif