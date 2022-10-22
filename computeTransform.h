#ifndef COMPUTETRANSFORM_H  
#define COMPUTETRANSFORM_H

#include <vector>
#include <iostream>
#include <string>
#include <chrono>
#include <unistd.h>
#include <stack>
#include <algorithm>

#include "lib/optSAIS.h" // algo for opt bwt construction 
#include "common.hpp"

#include "external/malloc_count/malloc_count.h"
#include "external/malloc_count/stack_count.h"

const int alph_size = 128;

void compute_optbwt_unsigned(Args arg);

void compute_bwt_sap_unsigned(Args arg); 

#endif