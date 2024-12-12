from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py as _build_py
import numpy
import os

os.environ['CC'] = 'g++'


class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()


compile_args = ['-std=c++17', '-O3', '-fPIC']  # , '-DNO_QC']
swig_args = ['-c++', '-builtin', '-doxygen', '-small', '-O', ]

constants_module = Extension('nwsspc.sharp.calc._constants',
                             sources=['nwsspc/sharp/calc/constants.i'],
                             swig_opts=swig_args,
                             extra_compile_args=compile_args
                             )

interp_module = Extension('nwsspc.sharp.calc._interp',
                          sources=['nwsspc/sharp/calc/interp.i',
                                   '../SHARPlib/interp.cpp'],
                          swig_opts=swig_args,
                          extra_compile_args=compile_args
                          )

# params_module = Extension('nwsspc.sharp.calc._params',
#         sources = ['nwsspc/sharp/calc/params.i', '../src/SHARPlib/parcel.cpp',
#                    '../src/SHARPlib/layer.cpp', '../src/SHARPlib/interp.cpp',
#                    '../src/SHARPlib/winds.cpp', '../src/SHARPlib/thermo.cpp',
#                    '../src/SHARPlib/params/convective.cpp'],
#         swig_opts = swig_args,
#         extra_compile_args = compile_args
#     )

parcel_module = Extension('nwsspc.sharp.calc._parcel',
                          sources=['nwsspc/sharp/calc/parcel.i',
                                   '../SHARPlib/layer.cpp',
                                   '../SHARPlib/interp.cpp',
                                   '../SHARPlib/thermo.cpp',
                                   '../SHARPlib/parcel.cpp'],
                          swig_opts=swig_args,
                          extra_compile_args=compile_args
                          )

profile_module = Extension('nwsspc.sharp.calc._profile',
                           sources=['nwsspc/sharp/calc/profile.i',
                                    '../SHARPlib/profile.cpp'],
                           swig_opts=swig_args,
                           extra_compile_args=compile_args
                           )

thermo_module = Extension('nwsspc.sharp.calc._thermo',
                          sources=['nwsspc/sharp/calc/thermo.i',
                                   '../SHARPlib/layer.cpp',
                                   '../SHARPlib/interp.cpp',
                                   '../SHARPlib/thermo.cpp'],
                          swig_opts=swig_args,
                          extra_compile_args=compile_args
                          )

winds_module = Extension('nwsspc.sharp.calc._winds',
                         sources=['nwsspc/sharp/calc/winds.i',
                                  '../SHARPlib/layer.cpp',
                                  '../SHARPlib/interp.cpp',
                                  '../SHARPlib/winds.cpp'],
                         swig_opts=swig_args,
                         extra_compile_args=compile_args
                         )

layer_module = Extension('nwsspc.sharp.calc._layer',
                         sources=['nwsspc/sharp/calc/layer.i',
                                  '../SHARPlib/interp.cpp',
                                  '../SHARPlib/layer.cpp'],
                         swig_opts=swig_args,
                         extra_compile_args=compile_args
                         )

setup(
    name="nwsspc_sharp_calc",
    version="0.0.1",
    cmdclass={'build_py': build_py},
    author="Kelton Halbert",
    author_email="kelton.halbert@noaa.gov",
    description="Used for processing and serving atmospheric sounding data.",
    url="https://github.com/keltonhalbert/SHARP-calc.git",
    classifiers=[
        "Programming Language :: Python :: 3",
                "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    packages=["nwsspc.sharp.calc"],
    include_dirs=[numpy.get_include(), "../../include",
                  "../../external/fmt/include"],
    ext_modules=[
        constants_module,
        interp_module,
        # params_module,
        parcel_module,
        profile_module,
        thermo_module,
        winds_module,
        layer_module
    ],
    py_modules=["nwsspc.sharp.calc.constants", "nwsspc.sharp.calc.interp",
                # "nwsspc.sharp.calc.params", "nwsspc.sharp.calc.parcel",
                "nwsspc.sharp.calc.parcel",
                "nwsspc.sharp.calc.profile", "nwsspc.sharp.calc.thermo",
                "nwsspc.sharp.calc.winds", "nwsspc.sharp.calc.layer"],
)
