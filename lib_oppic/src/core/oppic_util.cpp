#include <oppic_util.h>
#include <chrono>
#include <numeric>

//********************************************************************************
std::string getTimeStr()
{
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y_%m_%d__%H_%M_%S", std::localtime(&now));
    return s;
}

//********************************************************************************
std::vector<size_t> sort_indexes(const int* cell_indices, int size) 
{ 
    std::vector<size_t> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(), 
        [cell_indices](size_t i1, size_t i2) 
        { 
            return cell_indices[i1] < cell_indices[i2];
        });

    return idx;
}

//********************************************************************************



//********************************************************************************