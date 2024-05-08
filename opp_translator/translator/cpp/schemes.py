from pathlib import Path
from typing import Any, Dict

import cpp.translator.kernels as ctk
import op as OP
from language import Lang
from scheme import Scheme
from store import Application, ParseError, Program
from target import Target

class CppSeq(Scheme):
    lang = Lang.find("cpp")
    target = Target.find("seq")

    fallback = None

    consts_template = None
    loop_host_template = Path("cpp/seq/loop_host.hpp.jinja")
    move_loop_host_template = Path("cpp/seq/move_loop_host.hpp.jinja")
    master_kernel_template = Path("cpp/seq/master_kernel.cpp.jinja")
    header = "opp_seq.h"
    
    def translateKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        kernel_entities = app.findEntities(loop.kernel, program)
        if len(kernel_entities) == 0:
            raise ParseError(f"unable to find kernel: {loop.kernel}")

        extracted_entities = ctk.extractDependencies(kernel_entities, app)
        return ctk.writeSource(extracted_entities)

Scheme.register(CppSeq)

class CppMpi(Scheme):
    lang = Lang.find("cpp")
    target = Target.find("mpi")

    fallback = None

    consts_template = None
    loop_host_template = Path("cpp/mpi/loop_host.hpp.jinja")
    move_loop_host_template = Path("cpp/mpi/move_loop_host.hpp.jinja")
    master_kernel_template = Path("cpp/mpi/master_kernel.cpp.jinja")
    header = "opp_mpi.h"
    
    def translateKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        kernel_entities = app.findEntities(loop.kernel, program)
        if len(kernel_entities) == 0:
            raise ParseError(f"unable to find kernel: {loop.kernel}")

        extracted_entities = ctk.extractDependencies(kernel_entities, app)
        return ctk.writeSource(extracted_entities)

Scheme.register(CppMpi)

class CppOmp(Scheme):
    lang = Lang.find("cpp")
    target = Target.find("omp")

    fallback = None

    consts_template = None
    loop_host_template = Path("cpp/omp/loop_host.hpp.jinja")
    move_loop_host_template = Path("cpp/omp/move_loop_host.hpp.jinja")
    master_kernel_template = Path("cpp/omp/master_kernel.cpp.jinja")
    header = "opp_omp.h"
    
    def translateKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        kernel_entities = app.findEntities(loop.kernel, program)
        if len(kernel_entities) == 0:
            raise ParseError(f"unable to find kernel: {loop.kernel}")

        extracted_entities = ctk.extractDependencies(kernel_entities, app)
        return ctk.writeSource(extracted_entities)

# Scheme.register(CppOmp)

class CppCuda(Scheme):
    lang = Lang.find("cpp")
    target = Target.find("cuda")

    fallback = None

    consts_template = None
    loop_host_template = Path("cpp/cuda/loop_host.hpp.jinja")
    move_loop_host_template = Path("cpp/cuda/move_loop_host.hpp.jinja")
    master_kernel_template = Path("cpp/cuda/master_kernel.cu.jinja")
    header = "opp_cuda.h"
    
    def translateKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        kernel_entities = app.findEntities(loop.kernel, program)
        if len(kernel_entities) == 0:
            raise ParseError(f"unable to find kernel: {loop.kernel}")

        extracted_entities = ctk.extractDependencies(kernel_entities, app)

        if loop.loop_type == OP.LoopType.MOVE_LOOP:
            ctk.updateMoveKernelArgs(extracted_entities, 
                lambda typ, _: f"char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen\n    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen\n    {typ}", loop.kernel)

        ctk.updateFunctionTypes(extracted_entities, lambda typ, _: f"__device__ {typ}")
        ctk.renameConsts(extracted_entities, app, lambda const, _: f"{const}_d")

        for entity, rewriter in filter(lambda e: e[0] in kernel_entities, extracted_entities):
            kernel_idx = 0
            for i, (l, program) in enumerate(app.uniqueLoops(), 1):
                if l.kernel == loop.kernel:
                    kernel_idx = i
            ctk.insertStrides(
                entity,
                rewriter,
                app,
                loop,
                lambda dat_id: f"opp_k{kernel_idx}_{dat_id}_stride_d",
                skip=lambda arg: arg.access_type == OP.AccessType.INC and (arg.map_id is not None or arg.p2c_id is not None) and config["atomics"],
            )

        return ctk.writeSource(extracted_entities)

Scheme.register(CppCuda)

# class CppMPIOpenMP(Scheme):
#     lang = Lang.find("cpp")
#     target = Target.find("mpi_openmp")

#     const_template = None
#     loop_host_template = Path("cpp/mpi_openmp/loop_host.cpp.j2")
#     master_kernel_template = Path("cpp/mpi_openmp/master_kernel.cpp.j2")

#     loop_kernel_extension = "cpp"
#     master_kernel_extension = "cpp"

#     def translateKernel(
#         self, 
#         loop: ops.Loop, 
#         program: Program, 
#         app: Application, 
#         kernel_idx: int
#     ) -> str:
#         kernel_entities = app.findEntities(loop.kernel, program)

#         if len(kernel_entities) == 0:
#             raise ParseError(f"Unable to find kernel: {loop.kernel}")

#         extracted_entities = ctk.extractDependancies(kernel_entities, app)
#         return ctk.writeSource(extracted_entities)

# Scheme.register(CppMPIOpenMP)


# class CppCuda(Scheme):
#     lang = Lang.find("cpp")
#     target = Target.find("cuda")

#     const_template = None
#     loop_host_template = Path("cpp/cuda/loop_host.cpp.j2")
#     master_kernel_template = Path("cpp/cuda/master_kernel.cpp.j2")

#     loop_kernel_extension = "cu"
#     master_kernel_extension = "cu"

#     def translateKernel(
#         self,
#         loop: ops.Loop,
#         program: Program,
#         app: Application,
#         kernel_idx: int
#     ) -> str:
#         kernel_entities = app.findEntities(loop.kernel, program)

#         if len(kernel_entities) == 0:
#             raise ParseError(f"Unable to find kernel: {loop.kernel}")

#         extracted_entities = ctk.extractDependancies(kernel_entities, app)
#         return ctk.writeSource(extracted_entities)

# Scheme.register(CppCuda)


# class CppHip(Scheme):
#     lang = Lang.find("cpp")
#     target = Target.find("hip")

#     const_template = None
#     loop_host_template = Path("cpp/cuda/loop_host.cpp.j2")
#     master_kernel_template = Path("cpp/cuda/master_kernel.cpp.j2")

#     loop_kernel_extension = "cpp"
#     master_kernel_extension = "cpp"

#     def translateKernel(
#         self,
#         loop: ops.Loop,
#         program: Program,
#         app: Application,
#         kernel_idx: int
#     ) -> str:
#         kernel_entities = app.findEntities(loop.kernel, program)

#         if len(kernel_entities) == 0:
#             raise ParseError(f"Unable to find kernel: {loop.kernel}")

#         extracted_entities = ctk.extractDependancies(kernel_entities, app)
#         return ctk.writeSource(extracted_entities)

# Scheme.register(CppHip)


# class CppOpenMPOffload(Scheme):
#     lang = Lang.find("cpp")
#     target = Target.find("openmp_offload")

#     const_template = None
#     loop_host_template = Path("cpp/openmp_offload/loop_host.cpp.j2")
#     master_kernel_template = Path("cpp/openmp_offload/master_kernel.cpp.j2")

#     loop_kernel_extension = "cpp"
#     master_kernel_extension = "cpp"

#     def translateKernel(
#         self,
#         loop: ops.Loop,
#         program: Program,
#         app: Application,
#         kernel_idx: int
#     ) -> str:
#         kernel_entities = app.findEntities(loop.kernel, program)

#         if len(kernel_entities) == 0:
#             raise ParseError(f"Unable to find kernel: {loop.kernel}")

#         extracted_entities = ctk.extractDependancies(kernel_entities, app)
#         return ctk.writeSource(extracted_entities)

# Scheme.register(CppOpenMPOffload)


# #class CppOpenACC(Scheme):
# #    lang = Lang.find("cpp")
# #    target = Target.find("openacc")
# #
# #    const_template = None
# #    loop_host_template = Path("cpp/openacc/loop_host.cpp.j2")
# #    master_kernel_template = Path("cpp/openacc/master_kernel.cpp.j2")
# #
# #    loop_kernel_extension = "cpp"
# #    master_kernel_extension = "cpp"
# #
# #    def translateKernel(
# #        self,
# #        loop: ops.Loop,
# #        program: Program,
# #        app: Application,
# #        kernel_idx: int
# #    ) -> str:
# #        kernel_entities = app.findEntities(loop.kernel, program)
# #
# #        if len(kernel_entities) == 0:
# #            raise ParseError(f"Unable to find kernel: {loop.kernel}")
# #
# #        extracted_entities = ctk.extractDependancies(kernel_entities, app)
# #        return ctk.writeSource(extracted_entities)

# #Scheme.register(CppOpenACC)


# class CppSycl(Scheme):
#     lang = Lang.find("cpp")
#     target = Target.find("sycl")

#     const_template = None
#     loop_host_template = Path("cpp/sycl/loop_host.cpp.j2")
#     master_kernel_template = Path("cpp/sycl/master_kernel.cpp.j2")

#     loop_kernel_extension = "cpp"
#     master_kernel_extension = "cpp"

#     def translateKernel(
#         self,
#         loop: ops.Loop,
#         program: Program,
#         app: Application,
#         kernel_idx: int
#     ) -> str:
#         kernel_entities = app.findEntities(loop.kernel, program)

#         if len(kernel_entities) == 0:
#             raise ParseError(f"Unable to find kernel: {loop.kernel}")

#         extracted_entities = ctk.extractDependancies(kernel_entities, app)
#         return ctk.writeSource(extracted_entities)

# Scheme.register(CppSycl)
