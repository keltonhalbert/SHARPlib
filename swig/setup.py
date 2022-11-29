import setuptools

compile_args = ['-std=c++17']

interp_module = setuptools.Extension('nwsspc.sharp.calc._interp',
        sources = ['interp.i',
                   '../src/sharptab/profile.cpp', '../src/sharptab/utils.cpp',
                   '../src/sharptab/interp.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )


params_module = setuptools.Extension('nwsspc.sharp.calc._params',
        sources = ['params.i',
                   '../src/sharptab/profile.cpp', '../src/sharptab/utils.cpp',
                   '../src/sharptab/interp.cpp', '../src/sharptab/winds.cpp',
                   '../src/sharptab/thermo.cpp','../src/sharptab/params.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )


parcel_module = setuptools.Extension('nwsspc.sharp.calc._parcel',
        sources = ['parcel.i',
                   '../src/sharptab/profile.cpp', '../src/sharptab/utils.cpp',
                   '../src/sharptab/interp.cpp', '../src/sharptab/thermo.cpp',
                   '../src/sharptab/parcel.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )

profile_module = setuptools.Extension('nwsspc.sharp.calc._profile',
        sources = ['profile.i',
                   '../src/sharptab/interp.cpp', '../src/sharptab/utils.cpp',
                   '../src/sharptab/profile.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )

thermo_module = setuptools.Extension('nwsspc.sharp.calc._thermo',
        sources = ['thermo.i',
                   '../src/sharptab/utils.cpp',   '../src/sharptab/interp.cpp', 
                   '../src/sharptab/profile.cpp', '../src/sharptab/thermo.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )

winds_module = setuptools.Extension('nwsspc.sharp.calc._winds',
        sources=['winds.i',
                 '../src/sharptab/utils.cpp',   '../src/sharptab/interp.cpp', 
                 '../src/sharptab/profile.cpp', '../src/sharptab/winds.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )

utils_module = setuptools.Extension('nwsspc.sharp.calc._utils',
        sources = ['utils.i',
                   '../src/sharptab/profile.cpp', '../src/sharptab/interp.cpp',
                   '../src/sharptab/utils.cpp'],
        include_dirs=['../include'],
        swig_opts = ['-c++', '-builtin'],
        extra_compile_args = compile_args
    )

setuptools.setup(
	name="nwsspc_sharp_calc", 
	version="0.0.1",
	author="Kelton Halbert",
	author_email="kelton.halbert@noaa.gov",
	description="Used for processing and serving atmospheric sounding data.",
	url="https://github.com/keltonhalbert/SHARP-calc.git",
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: Apache 2.0 License",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.9',
    ext_modules = [interp_module, params_module, parcel_module, 
                   profile_module, thermo_module, winds_module, 
                   utils_module],
)
