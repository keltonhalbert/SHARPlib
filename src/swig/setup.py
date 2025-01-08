from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py as _build_py
import numpy
import os

os.environ['CC'] = 'g++'


class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()


compile_args = ['-std=c++17']
swig_args = ['-c++', '-builtin', '-doxygen', '-small', '-O', ]

constants_module = Extension(
    'nwsspc.sharp.calc._constants',
    sources=['nwsspc/sharp/calc/constants.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

interp_module = Extension(
    'nwsspc.sharp.calc._interp',
    sources=['nwsspc/sharp/calc/interp.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

params_module = Extension(
    'nwsspc.sharp.calc._params',
    sources=['nwsspc/sharp/calc/params.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

parcel_module = Extension(
    'nwsspc.sharp.calc._parcel',
    sources=['nwsspc/sharp/calc/parcel.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

thermo_module = Extension(
    'nwsspc.sharp.calc._thermo',
    sources=['nwsspc/sharp/calc/thermo.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

winds_module = Extension(
    'nwsspc.sharp.calc._winds',
    sources=['nwsspc/sharp/calc/winds.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
)

layer_module = Extension(
    'nwsspc.sharp.calc._layer',
    sources=['nwsspc/sharp/calc/layer.i'],
    swig_opts=swig_args,
    extra_objects=["../../build/libSHARPlib.a"]
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
        params_module,
        parcel_module,
        thermo_module,
        winds_module,
        layer_module
    ],
    py_modules=["nwsspc.sharp.calc.constants", "nwsspc.sharp.calc.interp",
                "nwsspc.sharp.calc.params", "nwsspc.sharp.calc.parcel",
                "nwsspc.sharp.calc.thermo", "nwsspc.sharp.calc.winds", 
                "nwsspc.sharp.calc.layer"],
)
