#include "main.h"
#include <string>

[[nodiscard]] void aggregate() noexcept
{
    size_t cnt_1 = std::stoi(input1All());
    size_t cnt_2 = std::stoi(input2All());
    output(std::to_string(cnt_1 + cnt_2));
}