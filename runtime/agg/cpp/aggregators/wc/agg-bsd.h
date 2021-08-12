#ifndef AGG_WC_H
#define AGG_WC_H

#include "common.h"

int platform_dependent_option_count_tweak(int x) { 
    return x - (g_options.is_present(2) && g_options.is_present(3)); 
}

#endif