#pragma once

#include "opp_types.h"

namespace opp::util
{
    template<typename T> 
    std::vector<T> getVecFromPtr(T* data, UINT size)
    {
        std::vector<T> dataVec(data, data + size);
        return dataVec;
    };  
};