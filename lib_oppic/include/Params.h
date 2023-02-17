#pragma once

#include <string>
#include <map>

#define INT int
#define REAL double
#define BOOL bool
#define STRING std::string

namespace opp
{
    class Params
    {
        public:
            Params(std::string fileanme);

            template <typename T>
            void add(std::string key, T value)
            {
                if constexpr (std::is_same<T, INT>::value)
                {
                    int_params.insert({key, value});
                }
                else if constexpr (std::is_same<T, REAL>::value)
                {
                    real_params.insert({key, value});
                }
                else if constexpr (std::is_same<T, BOOL>::value)
                {
                    bool_params.insert({key, value});
                }
                else if constexpr (std::is_same<T, STRING>::value)
                {
                    str_params.insert({key, value});
                }

            }

            template <typename T>
            T get(std::string key)
            {
                if constexpr (std::is_same<T, INT>::value)
                {
                    return int_params[key];
                }
                else if constexpr (std::is_same<T, REAL>::value)
                {
                    return real_params[key];
                }
                else if constexpr (std::is_same<T, BOOL>::value)
                {
                    return bool_params[key];
                }
                else if constexpr (std::is_same<T, STRING>::value)
                {
                    return str_params[key];
                }

                return NULL;
            }

            void write(std::ostream &out);

        private:
            std::map<std::string, STRING> str_params;
            std::map<std::string, INT> int_params;
            std::map<std::string, REAL> real_params;
            std::map<std::string, BOOL> bool_params;

    };
    
};
