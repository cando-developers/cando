/*
    File: foundation.h
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
#ifndef cando_foundation_H
#define cando_foundation_H

#define Chem_temp_Defun(x) core::af_def(ChemPkg, #x, &chem__##x, ARGS_chem__##x, DECL_chem__##x, DOCS_chem__##x, __FILE__, __LINE__);
#define Adapt_temp_Defun(x) core::af_def(AdaptPkg, #x, &adapt__##x, ARGS_adapt__##x, DECL_adapt__##x, DOCS_adapt__##x, __FILE__, __LINE__);
#define Geom_temp_Defun(x) core::af_def(GeomPkg, #x, &geom__##x, ARGS_geom__##x, DECL_geom__##x, DOCS_geom__##x, __FILE__, __LINE__);
#define Units_temp_Defun(x) core::af_def(UnitsPkg, #x, &units__##x, ARGS_units__##x, DECL_units__##x, DOCS_units__##x, __FILE__, __LINE__);

#endif
