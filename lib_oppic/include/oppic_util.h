#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <cstring>

//********************************************************************************
std::string getTimeStr();

//********************************************************************************
inline char *copy_str(char const *src) 
{
    char *dest = (char *)malloc((strlen(src) + 1) * sizeof(char));
    return strncpy(dest, src, strlen(dest));
}

//********************************************************************************
std::vector<size_t> sort_indexes(const int* cell_indices, int size);

//********************************************************************************