from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(
        ["crhs.pyx", "numerov.pyx", "poisson.pyx", "numerovgen.pyx", "argums.pyx", "legendre.pyx", 
        "ho.pyx", "root_mu.pyx", "compute_mt.pyx"], 
        language_level="3"
        ),
    include_dirs=[numpy.get_include()]
)
