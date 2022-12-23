#pragma once

#include <oppic_lib.h>


template <typename... T, typename... OPARG>
void oppic_par_loop(void (*kernel)(T *...), char const *name, oppic_set set, oppic_iterate_type iter_type,
                 OPARG... arguments) {
    printf("oppic_par_loop %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}

template <typename... T, typename... OPARG>
void oppic_par_loop_particle(void (*kernel)(T *...), char const *name, oppic_set set, oppic_iterate_type iter_type,
                 OPARG... arguments) {
    printf("oppic_par_looppic_particle %s iterate %s\n", name, (iter_type == OP_ITERATE_ALL) ? "all" : "only injected");
}

