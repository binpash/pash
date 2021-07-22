#ifndef MAIN_H
#define MAIN_H

// See instructions below this code
#include <string>
#include <istream>
#include "command_opts.h"

// -*-*-*-*-*-*-*-*-*-* INSTRUCTIONS BEGIN HERE -*-*-*-*-*-*-*-*-*-*
// This file contains common code between all aggregators
// It should be included when implementing a new aggregator
// It contains the aggregator's main function
// To use it, just include it and implement the function presented below
// and remember to link the executable with main.cpp

void aggregate() noexcept; // implement this

// It should work in the following manner:
// 1. The aggregator should create a global variable called
//    g_options in the manner showed below:

// constexpr cmd_opts g_options{
//     cmd_opt{ "combine", 'c'},
//     cmd_opt{ "reverse", 'r'}
// };

//   Then, the member functions is_present and argument can be used
//   with the arguemnt being the index of the option in the constructor.

// 2. It should request inputs to aggregate when it needs them using
//    objects returned by the functions below:

[[nodiscard]] std::istream& input1() noexcept;
[[nodiscard]] std::istream& input2() noexcept;

// 3. It should output the aggregated results as soon as it has a piece
//    of them done to the stream returned to function below:

[[nodiscard]] std::ostream& output() noexcept;

// 4. If some behaviour has not been implemented yet, the below
//    function can be called:

[[noreturn]] void nyi_error(const char* message) noexcept;


#endif