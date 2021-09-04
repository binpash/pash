#ifndef AGG_WC_H
#define AGG_WC_H

#include "common.h"

size_t platform_dependent(int&) {
    size_t padding = std::numeric_limits<size_t>::max();
    {
        size_t prev_pos = 0, dummy;
        while(input1() >> dummy)
        {
            size_t new_pos = input1().tellg();
            padding = std::min(new_pos - prev_pos, padding);
            prev_pos = new_pos + 1;
        }
        input1().clear();
        input1().seekg(std::ios_base::beg);
    }
    return padding;
}

#endif
