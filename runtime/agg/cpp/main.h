// See instructions below this include block
#include <string>
#include <optional>
#include <string_view>
#include <vector>

// -*-*-*-*-*-*-*-*-*-* INSTRUCTIONS BEGIN HERE -*-*-*-*-*-*-*-*-*-*
// This file contains common code between all aggregators
// It should be included when implementing a new aggregator
// It contains the aggregator's main function
// To use it, just include it and implement the function presented below
// and remember to link the executable with main.cpp

[[nodiscard]] void aggregate() noexcept; // implement this

// It should work in the following manner:
// 1. It should request inputs to aggregate when it needs them using the functions below:

// Read a single line from input. Returns std::nullopt if input has been exhausted
[[nodiscard]] std::optional<std::string> input1() noexcept;
[[nodiscard]] std::optional<std::string> input2() noexcept;

// Read the entire input as one large string, possibly containing multiple lines
[[nodiscard]] std::string input1All() noexcept;
[[nodiscard]] std::string input2All() noexcept;

// Read the entire input as one large string. Returns each line in a separate string_view.
// The first element of the pair has to be saved as it owns the memory pointed to by the string_views
[[nodiscard]] std::pair<std::string, std::vector<std::string_view>> input1AllLines() noexcept;
[[nodiscard]] std::pair<std::string, std::vector<std::string_view>> input2AllLines() noexcept;

//    These functions block until input is available in the appropriate file and
//    return a single line as soon as it is available.
// 2. It should call to function below to output the aggregated results as soon as it has a piece
//    of it done. It can either be a single line, multiple lines or the entire output at once.
//    This function appends a '\n' if the given string doesn't end in it

void output(const std::string&) noexcept;
void output(const char*, size_t) noexcept;
