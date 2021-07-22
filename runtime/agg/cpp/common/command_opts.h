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

using cmd_opt = std::pair<const char*, char>;

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
                    if (opt.second)
                    {
                        map[opt.second] = i;
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
                if (opt.first)
                {
                    long_opts[current_idx++] = {
                        opt.first,
                        optional_argument,
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
                if (opt.second)
                {
                    optstring[current_idx++] = opt.second;
                    optstring[current_idx++] = ':';
                    optstring[current_idx++] = ':';
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