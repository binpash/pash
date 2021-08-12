#ifndef AGG_WC_H
#define AGG_WC_H

#include "common.h"

size_t platform_dependent(int& numbers_to_input) { 
    if(g_options.is_present(2) && g_options.is_present(3))
        --numbers_to_input;
    return 7;
}

#endif
