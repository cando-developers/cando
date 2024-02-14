/*
    File: extension.fwd.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#ifndef extension_fwd_H
#define extension_fwd_H

#include <cstddef>
#include <cstdint>

#define VECREAL_DOUBLE 2
#define VECREAL_FLOAT  1

//
// This is where you define the VECREAL and NUM_REAL types
//
// FIXME: Switch back to VECREAL VECREAL_FLOAT
#define VECREAL VECREAL_DOUBLE
#define NUM_REAL VECREAL_DOUBLE

#if VECREAL == VECREAL_DOUBLE
typedef double  vecreal;
typedef int64_t vecint;
#define NVECTOR_O_TYPE core::SimpleVector_double_O
#define NVECTOR_sp_TYPE core::SimpleVector_double_sp
#else
typedef float vecreal;
typedef int32_t vecint;
#define NVECTOR_O_TYPE core::SimpleVector_float_O
#define NVECTOR_sp_TYPE core::SimpleVector_float_sp
#endif

#if NUM_REAL == VECREAL_DOUBLE
typedef double num_real;
typedef int64_t num_int;
#else
typedef float  num_real;
typedef int32_t num_int;
#endif

extern double safe_acos(double rad);
extern double safe_asin(double rad);
extern float  safe_acos(float rad);
extern float  safe_asin(float rad);


namespace cando {
// Store the number of doubles that fit into a SIMD register
// global_simd_width
//     1 = No SIMD
//     2 doubles SSE
//     4 doubles AVX2
//     8 doubles AVX-512
//
extern size_t global_simd_width;

};
#endif
