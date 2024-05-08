import re

from store import Program
from util import SourceBuffer
import op as OP

# Augment source program to use generated kernel hosts
def translateProgram(source: str, program: Program, force_soa: bool) -> str:
    buffer = SourceBuffer(source)

    # 1. Remove the defined constants : OPP_* CONST_*
    pattern = re.compile(r'^\s*(OPP_REAL|OPP_BOOL|OPP_INT)\sCONST_.*$')
    lines_removed = 0
    index = buffer.search(pattern)
    while index is not None:
        buffer.remove(index)
        source = buffer.translate()
        buffer = SourceBuffer(source)
        lines_removed += 1
        index = buffer.search(pattern)

    # 2. Update loop calls
    for loop in program.loops:

        print(f'loop {loop.kernel}')

        if loop.loop_type == OP.LoopType.PAR_LOOP:
            before, after = buffer.get(loop.loc.line - lines_removed - 1).split("opp_par_loop", 1)
            after = re.sub(
                rf"{loop.kernel}\s*,\s*\"{loop.loop_name}\"\s*,\s*", "", after, count=1
            )  # TODO: This assumes that the kernel arg is on the same line as the call
            buffer.update(loop.loc.line - lines_removed - 1, before + f"opp_par_loop_{loop.iterator_type.name}__{loop.kernel}" + after)

        elif loop.loop_type == OP.LoopType.MOVE_LOOP:
            before, after = buffer.get(loop.loc.line - lines_removed - 1).split("opp_particle_move", 1)
            after = re.sub(
                rf"{loop.kernel}\s*,\s*\"{loop.loop_name}\"\s*,\s*", "", after, count=1
            )  # TODO: This assumes that the kernel arg is on the same line as the call
            buffer.update(loop.loc.line - lines_removed - 1, before + f"opp_particle_move__{loop.kernel}" + after)

    # 3. Update headers
    index = buffer.search(r'\s*#include\s+"opp_templates\.h"')
    assert index is not None
    index += 2

    # buffer.insert(index, '#ifdef OPENACC\n#ifdef __cplusplus\nextern "C" {\n#endif\n#endif\n')
    loop_list = []
    for loop in program.loops:
        if loop.loop_type == OP.LoopType.PAR_LOOP:
            prototype = f'void opp_par_loop_{loop.iterator_type.name}__{loop.kernel}(opp_set,opp_iterate_type{",opp_arg" * len(loop.args)});'
        elif loop.loop_type == OP.LoopType.MOVE_LOOP:
            prototype = f'void opp_particle_move__{loop.kernel}(opp_set,opp_map,opp_dat{",opp_arg" * len(loop.args)});'
        
        if prototype not in loop_list:
            loop_list.append(prototype)
            buffer.insert(index, prototype)

    source = buffer.translate()

    # Substitude the opp_templates.h include for opp_lib.h
    source = re.sub(r'#include\s+"opp_templates\.h"', '#include "opp_lib.h"', source)
    source = re.sub(r'// USER WRITTEN CODE', '// AUTO GENERATED CODE', source)
    source = re.sub(r'#include "kernels.h"', '// #include "kenels.h" // codegen commented :  TODO: ...', source)

    return source
