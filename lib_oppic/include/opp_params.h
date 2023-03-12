#pragma once

#include <string>
#include <map>
#include <limits.h>
#include <float.h>

#ifndef INT
    #define INT int
    #define REAL double
    #define BOOL bool
    #define STRING std::string
#endif

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace opp
{
    class Params
    {
        public:
            Params(std::string fileanme);

            template <typename T>
            inline void add(std::string key, T value) // refer specilized template functions
            {
                UNUSED(value);
                std::cerr << "Params add not implemented for type of : " << key << std::endl;
                // exit(-1);
            }

            template <typename T>
            inline T get(std::string key) // refer specilized template functions
            {
                std::cerr << "Params get not implemented for type of : " << key << std::endl;
                // exit(-1);
                return NULL;
            }

            void write(std::ostream &out);

        private:
            std::map<std::string, STRING> str_params;
            std::map<std::string, INT> int_params;
            std::map<std::string, REAL> real_params;
            std::map<std::string, BOOL> bool_params;
    };
    
    template <>
    inline void Params::add<INT>(std::string key, INT value)
    {
        int_params[key] = value;
    }

    template <>
    inline void Params::add<REAL>(std::string key, REAL value)
    {
        real_params[key] = value;
    }

    template <>
    inline void Params::add<BOOL>(std::string key, BOOL value)
    {
        bool_params[key] = value;
    }

    template <>
    inline void Params::add<STRING>(std::string key, STRING value)
    {
        str_params[key] = value;
    }

    template <>
    inline INT Params::get<INT>(std::string key)
    {
        auto a = int_params.find(key);
        if (a != int_params.end()) return a->second;
        
        std::cerr << "ERROR : " << key << " (INT) not found in configs" << std::endl;
        return INT_MAX; 
    }

    template <>
    inline REAL Params::get<REAL>(std::string key)
    {
        auto a = real_params.find(key);
        if (a != real_params.end()) return a->second;
        
        std::cerr << "ERROR : " << key << " (REAL) not found in configs" << std::endl;
        return DBL_MAX; 
    }
    
    template <>
    inline BOOL Params::get<BOOL>(std::string key)
    {
        auto a = bool_params.find(key);
        if (a != bool_params.end()) return a->second;
        
        std::cerr << "ERROR : " << key << " (BOOL) not found in configs" << std::endl;
        return false; 
    }

    template <>
    inline STRING Params::get<STRING>(std::string key)
    {
        auto a = str_params.find(key);
        if (a != str_params.end()) return a->second;
        
        std::cerr << "ERROR : " << key << " (STRING) not found in configs" << std::endl;
        return "NOT FOUND"; 
    }
};
