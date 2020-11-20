import os
from distutils.core import setup, Extension
from Cython.Build import cythonize

mpi_compile_args = os.popen("mpic++ --showme:compile").read().strip().split(' ')
mpi_link_args    = os.popen("mpic++ --showme:link").read().strip().split(' ')

mpi_compile_args.extend(['-O3','-fopenmp'])
mpi_link_args.extend(['-lgomp',])
os.environ['CFLAGS'] = ' '.join(mpi_compile_args)
os.environ['LDFLAGS'] = ' '.join(mpi_link_args)

ext = Extension('epiclib',
                sources=["epiclib.pyx","../../../EPIC/beambeam.cpp", "../../../EPIC/acc_model.cpp", "../../../EPIC/beam.cpp",
                    "../../../EPIC/faddeeva.cpp", "../../../EPIC/mathfunc.cpp",],  # additional source file(s)

                extra_compile_args = mpi_compile_args.extend(['-O3','-fopenmp']),
                extra_link_args    = mpi_link_args,
                language="c++",
                )

setup (name='epiclib', ext_modules=cythonize(ext))

