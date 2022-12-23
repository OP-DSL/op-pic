/*==============================================================================*
 * TRACE
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-08-24
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
// #include <omp.h>

#define USE_TRACE

#ifdef USE_TRACE


#define TRACE_ME TraceCaller _TRACE_OBJECT(__func__);
#define TRACE(x) TraceCaller _TRACE_OBJECT(x);

double _walltime();

struct TraceElement {
    std::string name;
    int num_calls = 0;
    double total_time = 0;
    double child_time = 0;
    double start_time = 0;
};

class Trace {
    private:
        double entry_time = 0;
        double exit_time = 0;
        double trace_start;
        std::string name;
        std::vector<TraceElement*> callstack;
        std::unordered_map<std::string, TraceElement> function_list;

    public:
        Trace(std::string name_);

        void write_callstack();
        void write_profile(std::string filename);
        void enter(std::string func_name);
        void exit(std::string func_name);
};

namespace trace {
    extern int enabled;
    extern Trace current;
}

#else /* !USE_TRACE */

    #define TRACE_ME
    #define TRACE(x)
#endif 

struct TraceCaller {
    std::string name;
    TraceCaller(std::string name_);
    ~TraceCaller();

    void add(std::string name_);
};

#endif /* !TRACE_H */
