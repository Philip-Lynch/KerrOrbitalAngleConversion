// common_definitions.h
#ifndef COMMON_DEFINITIONS_H
#define COMMON_DEFINITIONS_H

#include <math.h>
// Definitions for using CForm'ed output from Mathematica
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Pi				M_PI

#endif // COMMON_DEFINITIONS_H