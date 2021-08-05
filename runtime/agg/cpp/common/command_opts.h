#ifndef COMMAND_OPTS_H
#define COMMAND_OPTS_H

#include <vector>
#include <array>
#include <string>
#include <initializer_list>
#include <getopt.h>

struct opt_holder
{
    option* long_options;
    const char* optstring;
    int* present;
    std::string* args;
    std::array<char, 256> map;
};

struct cmd_opt
{
    const char* long_name;
    char abbreviation;
    enum class Argument {
        None = no_argument,
        Required = required_argument,
        Optional = optional_argument
    } argument;
};

template <size_t Size>
class cmd_opts : public std::array<cmd_opt, Size>, public opt_holder
{
public:
    constexpr cmd_opts(std::initializer_list<cmd_opt> ar) noexcept :
        std::array<cmd_opt, Size>{il_to_array<Size>(ar)},
        opt_holder{long_opts.data(), optstring.data(), present, arg,
            [&](){
                std::array<char, 256> map{};
                auto& base = static_cast<std::array<cmd_opt, Size>&>(*this);
                int i = 0;
                for (auto& opt : base)
                {
                    if (opt.abbreviation)
                    {
                        map[opt.abbreviation] = i;
                    }
                    ++i;
                }
                return map;
            }()
        },
        long_opts{[&](){
            std::array<option, Size + 1> long_opts{};
            size_t current_idx = 0;
            auto& base = static_cast<std::array<cmd_opt, Size>&>(*this);
            int i = 0;
            for (auto& opt : base)
            {
                if (opt.long_name && opt.long_name[0] != '\0')
                {
                    long_opts[current_idx++] = {
                        opt.long_name,
                        static_cast<int>(opt.argument),
                        &present[i],
                        1
                    };
                }
                ++i;
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
                    if (opt.argument != cmd_opt::Argument::None)
                    {
                        optstring[current_idx++] = ':';
                        if (opt.argument == cmd_opt::Argument::Optional)
                        {
                            optstring[current_idx++] = ':';
                        }
                    }
                }
            }
            return optstring;
        }()}
        {
        }

        [[nodiscard]] bool is_present(size_t opt_idx) const noexcept
        {
            return present[opt_idx] == 1;
        }
        [[nodiscard]] const std::string& argument(size_t opt_idx) const noexcept
        {
            return arg[opt_idx];
        }
private:
    std::array<option, Size + 1> long_opts;
    std::array<char, 3 * Size + 1> optstring;
    inline static int present[Size + 1];
    inline static std::string arg[Size + 1];

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
