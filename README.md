# KerrOrbitalAngleConversion
A repository for converting between the different orbital angles used to describe bound timelike geodesics in Kerr spacetime. This code contains a Wolfram Language (Mathematica) implimentation, a C implimentation and Cython bindings so that one can use the C implimentation in Python. 

## Mathematica Implimentation
The mathematica implimentation is located in the mma directory. Simply run any cell in the notebook and agree to evaluate all initialization cells. From there one can work through the "Testing the Conversion" section. 

## C Implimentation
This code requires the gsl library. 
To compile the C code, simply run the following command in the terminal:
```console
make
```

 Note that the Makefile has been configured to work with my M1 Macbook and may need to be changed to work on other machines.

This will create a binary file in the bin/ directory called "KerrAngleConversion" which takes 7 arguments in the following order: a p e x psi chi phi

Example:
```console
./KerrAngleConversion 0.9 10 0.3 0.8 1.2 0.5 10
```

## Python Implimentation
To create the python bindings using Cython, run the follwing command in the terminal: 
```console
python setup.py build_ext
```

To use, follow the example given in examples/example.ipynb

