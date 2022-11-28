from distutils.core import setup, Extension

compile_args = ['std-c++17']
swig_files = ['thermo.i', 'winds.i']

thermo_module = Extension('_sharpcalc_thermo',
        sources = ['thermo.i',
                   '../sharptab/utils.cpp',   '../sharptab/interp.cpp', 
                   '../sharptab/profile.cpp', '../sharptab/thermo.cpp'],
        include_dirs=['../../include'],
        swig_opts = ['-c++'],
        extra_args = compile_args
    )

winds_module = Extension('_sharpcalc_winds',
        sources=['winds.i',
                 '../sharptab/utils.cpp',   '../sharptab/interp.cpp', 
                 '../sharptab/profile.cpp', '../sharptab/winds.cpp'],
        include_dirs=['../../include'],
        swig_opts = ['-c++'],
        extra_args = compile_args
)

setup (name = 'nwsspc.sharp.calc',
   version = '0.1',
   author = 'Kelton Halbert',
   description = """Simple swig example from docs""",
   ext_modules = [thermo_module, winds_module],
   py_modules = ['nwsspc.sharp.calc.thermo', 'nwsspc.sharp.calc.winds']
)
