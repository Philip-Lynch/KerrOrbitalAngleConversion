from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# Define the extension module
ext_modules = [
    Extension(
        "kerrangleconversions",                               # Name of the module
        sources=[                                # Cython and C source files
            "cython/angle_conversions.pyx",                       # The Cython wrapper
            "src/angle_conversions.c",                     # Main C file
            "src/constants_of_motion.c",
            "src/oscillating_terms.c",
            "src/root_finding.c" 
        ],
        include_dirs=["include", "/usr/local/include","/opt/homebrew/opt/gsl/include", numpy.get_include()],  # Include directories
        library_dirs=["/usr/local/lib","/opt/homebrew/opt/gsl/lib"],                # Library directories
        libraries=["m", "gsl", "gslcblas"],                       # Link against libm (math) and libmylib
        extra_compile_args=["-O3"],                     # Compiler optimization
        extra_link_args=[]                              # Extra link arguments if needed
    )
    ]

# Build the Cython extension
setup(
    name="kerrangleconversions",
    ext_modules=cythonize(ext_modules),
)