#include <fstream>
#include <iostream>
#include <string>
#include <regex>

#include "Params.h"

using namespace opp;

Params::Params(std::string file_name) 
{
    std::ifstream params_file;
    params_file.open(file_name);

    std::string line;
    while(std::getline(params_file, line)) 
    {
        std::regex_replace(line, std::regex("\\s*\\(#.*\\)\\?$"), "");
        if (line.length() == 0) 
            continue;

        std::smatch matches;
        if (std::regex_match(line, matches, std::regex("(BOOL|STRING|INT|REAL)(\\s+)([aA-zZ]*)(\\s+)(\\=)(\\s+)(.*)"))) 
        {
            if (matches.size() != 8) 
            {
                std::cerr << "Invalid input line: '" << line << "'" << matches.size() << std::endl;
                exit(-1);
            }

            std::string type = matches[1];
            std::string key = matches[3];
            std::string value = matches[7];

            if (type == "STRING")
            {
                add<STRING>(key, value);
            }
            else if (type == "INT")
            {
                add<INT>(key, std::stoi(value));
            }
            else if (type == "REAL")
            {
                add<REAL>(key, std::stod(value));
            }
            else if (type == "BOOL")
            {
                bool insertValue = false;
                if (std::regex_match(value, std::regex("^t(rue)?$", std::regex_constants::icase)) || value == "1") 
                {
                    insertValue = true;
                } else if (std::regex_match(value, std::regex("^f(alse)?$", std::regex_constants::icase)) || value == "0") 
                {
                    insertValue = false;
                }

                add<BOOL>(key, insertValue);
            }
            else
            {
                std::cerr << "Unrecognised parameter: '" << line << "'" << std::endl;
                // exit(-1);
            }
        }

    }
}

void Params::write(std::ostream &out) {

    out << "SIMULATION PARAMETERS"  << std::endl;
    out << "---------------------" << std::endl;

    for (auto &a : str_params) {
        out << "    " << a.first << " = " << a.second << " (STRING)" << std::endl;
    }

    for (auto &a : int_params) {
        out << "    " << a.first << " = " << a.second << " (INT)" << std::endl;
    }

    for (auto &a : real_params) {
        out << "    " << a.first << " = " << a.second << " (REAL)" << std::endl;
    }

    for (auto &a : bool_params) {
        out << "    " << a.first << " = " << a.second << " (BOOL)" << std::endl;
    }

    out << "---------------------" << std::endl << std::endl << std::endl;
}

// int main()
// { Tests
//     Params p("/ext-home/zl/phd/test_any/system.param");
    
//     std::cout << "inlet_mesh : " << p.get<STRING>("inlet_mesh") << std::endl;
//     std::cout << "dt : " << p.get<REAL>("dt") << std::endl;
//     std::cout << "max_iter : " << p.get<INT>("max_iter") << std::endl;

//     p.write(std::cout);
// }