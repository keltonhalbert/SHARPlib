import setuptools

compile_args = ['-std=c++17']

thermo_module = setuptools.Extension('nwsspc.sharp.calc._thermo',
        sources = ['thermo.i',
                   '../sharptab/utils.cpp',   '../sharptab/interp.cpp', 
                   '../sharptab/profile.cpp', '../sharptab/thermo.cpp'],
        include_dirs=['../../include'],
        swig_opts = ['-c++'],
        extra_compile_args = compile_args
    )

winds_module = setuptools.Extension('nwsspc.sharp.calc._winds',
        sources=['winds.i',
                 '../sharptab/utils.cpp',   '../sharptab/interp.cpp', 
                 '../sharptab/profile.cpp', '../sharptab/winds.cpp'],
        include_dirs=['../../include'],
        swig_opts = ['-c++'],
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
    ext_modules = [thermo_module, winds_module],
)
