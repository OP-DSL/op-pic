from __future__ import annotations

import sys
import traceback
import re
from abc import abstractmethod
from pathlib import Path 
from typing import Any, Dict, List, Optional, Set, Tuple

from jinja2 import Environment

import op as OP
# import ops
from language import Lang
from store import Application, Program
from target import Target
from util import sycl_set_flat_parallel
from util import extract_intrinsic_functions
from util import Findable
from util import KernelProcess

class Scheme(Findable):
    lang: Lang
    target: Target

    loop_host_template: Path
    move_loop_host_template: Path
    master_kernel_template: Optional[Path]

    def __str__(self) -> str:
        return f"{self.lang.name}/{self.target.name}"

    def canGenLoopHost(self, loop: OP.Loop) -> bool:
        return True

    def getBaseConfig(self, loop: OP.Loop) -> Dict[str, Any]:
        return self.target.defaultConfig()

    def getConfig(self, loop: OP.Loop, config_overrides: List[Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
        config = self.getBaseConfig(loop)

        for override in config_overrides:
            for loop_match, items in override.items():
                if re.match(loop_match, loop.name):
                    config |= items

        return config
        
    def genLoopHost(
        self,
        env: Environment,
        loop: OP.Loop,
        program: Program,
        app: Application,
        kernel_idx: int,
        config_overrides: List[Dict[str, Dict[str, Any]]],
        force_generate: bool = False,
    ) -> Optional[Tuple[str, str, bool]]:
        
        if loop.loop_type == OP.LoopType.PAR_LOOP:
            template = env.get_template(str(self.loop_host_template))
            extension = self.loop_host_template.suffixes[-2][1:]
        else:
            template = env.get_template(str(self.move_loop_host_template))
            extension = self.move_loop_host_template.suffixes[-2][1:]            
        # Load the loop host template

        args = {
            "OP": OP,
            "lh": loop,
            "kernel_idx": kernel_idx,
            "lang": self.lang,
            "config": self.getConfig(loop, config_overrides),
            "kernel_func": None,
            "host_kernel_func": None,
        }

        try:
            if (not loop.fallback and self.canGenLoopHost(loop)) or force_generate:
                args["host_kernel_func"] = self.translateHostKernel(loop, program, app, args["config"], kernel_idx)
                args["kernel_func"] = self.translateKernel(loop, program, app, args["config"], kernel_idx)

        except Exception as e:
            print(f"Error: kernel translation for kernel {kernel_idx} ({loop.name}) failed ({self}):")
            print(f"  fallback: {loop.fallback}, can generate: {self.canGenLoopHost(loop)}, force_generate: {force_generate}")
            traceback.print_exc(file=sys.stdout)
        
        # print(f'ZAM genLoopHost AAAAAAAAAAA')
        # print(f'ZAM genLoopHost AAAAAAAAAAA {str(self.loop_host_template)}')


        if args["kernel_func"] is None and self.fallback is None:
            # print('if args["kernel_func"] is None and self.fallback is None')
            return None

        if self.fallback is not None:
            fallback_wrapper_template = env.get_template(str(self.lang.fallback_wrapper_template))
            fallback_template = env.get_template(str(self.fallback.loop_host_template))

            fallback_args = dict(args)

            fallback_args["config"] = self.fallback.getConfig(loop, config_overrides)
            fallback_args["kernel_func"] = self.fallback.translateKernel(
                loop, program, app, fallback_args["config"], kernel_idx
            )

        if args["kernel_func"] is None:
            # print('if args["kernel_func"] is None')
            return (fallback_template.render(**fallback_args, variant=""), extension, True)

        if self.fallback is None:
            # print('if self.fallback is None:')
            return (template.render(**args, variant=""), extension, False)

        source = template.render(**args, variant="_main")

        source += "\n\n"
        source += fallback_template.render(**fallback_args, variant="_fallback")

        source += "\n\n"
        source += fallback_wrapper_template.render(**args)

        return (source, extension, False)


        # template = env.get_template(str(self.loop_host_template))
        # #extention = self.loop_host_template.suffixes[-2][1:]

        # print(f'ZAM genLoopHost AAAAAAAAAAA')

        # kernel_func = self.translateKernel(loop, program, app, kernel_idx)

        # kp_obj = KernelProcess()
        # if(self.lang.name == "C++"):
        #     kernel_func = kp_obj.clean_kernel_func_text(kernel_func)
        
        #     if(self.target.name == "cuda"):
        #         kernel_func = kp_obj.cuda_complex_numbers(kernel_func)

        #     consts_in_kernel = []
        #     if(self.target.name == "sycl"):
        #         kernel_func, consts_in_kernel = kp_obj.sycl_kernel_func_text(kernel_func, app.consts())

        #     #TODO : Complex arguments in HIP

        #     const_dims = []
        #     if(self.target.name == "openacc"):
        #         consts_in_kernel, const_dims = kp_obj.openacc_get_const_names_and_dim(kernel_func, app.consts())

        #     kernel_body, args_list = kp_obj.get_kernel_body_and_arg_list(kernel_func)
        #     flat_parallel, ops_cpu = sycl_set_flat_parallel(loop.has_reduction)
        #     intrinsic_funcs = ""

        # elif (self.lang.name == "Fortran"):
        #     kernel_body = None
        #     consts_in_kernel = None
        #     const_dims = None
        #     args_list = None
        #     flat_parallel = None
        #     ops_cpu = None
        #     intrinsic_funcs = extract_intrinsic_functions(kernel_func)

        # # Generalte source from the template
        # return (
        #     template.render (
        #         ops=ops,
        #         lh=loop,
        #         kernel_func=kernel_func,
        #         kernel_idx=kernel_idx,
        #         kernel_body=kernel_body,
        #         consts_in_kernel=consts_in_kernel,
        #         const_dims=const_dims,
        #         args_list=args_list,
        #         intrinsic_funcs=intrinsic_funcs,
        #         lang=self.lang,
        #         target=self.target,
        #         soa_set=force_soa,
        #         flat_parallel=flat_parallel,
        #         ops_cpu=ops_cpu
        #     ),
        #     self.loop_kernel_extension
        # )

    def genMasterKernel(self, env: Environment, app: Application, user_types_file: Optional[Path]) -> Tuple[str, str]:
        if self.master_kernel_template is None:
            exit(f"No master kernel template registered for {self}")

        user_types = None
        if user_types_file is not None:
            user_types = user_types_file.read_text()

        # Load the loop host template
        template = env.get_template(str(self.master_kernel_template))

        extension = self.master_kernel_template.suffixes[-2][1:]
        name = f"opp_kernels.{extension}"

        unique_loops = []
        for loop, _ in app.loops():
            if loop.kernel not in unique_loops:
                unique_loops.append(loop.kernel)

        # print(f'at GEN MASTER KERNEL {app.consts()}')
        # Generate source from the template
        return template.render(OP=OP, app=app, lang=self.lang, target=self.target, user_types=user_types, header=self.header, unique_loops=unique_loops), name

        
    @abstractmethod
    def translateKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        pass
    
    def translateHostKernel(
        self,
        loop: OP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        return "// Host kernel not used"

    def matches(self, key: Tuple[Lang, Target]) -> bool:
        if not (isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], Lang) and isinstance(key[1], Target)):
            return False

        return self.lang == key[0] and self.target == key[1]