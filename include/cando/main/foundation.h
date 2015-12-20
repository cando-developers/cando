#ifndef cando_foundation_H
#define cando_foundation_H

#define Chem_temp_Defun(x) core::af_def(ChemPkg, #x, &chem__##x, ARGS_chem__##x, DECL_chem__##x, DOCS_chem__##x, __FILE__, __LINE__);
#define Adapt_temp_Defun(x) core::af_def(AdaptPkg, #x, &adapt__##x, ARGS_adapt__##x, DECL_adapt__##x, DOCS_adapt__##x, __FILE__, __LINE__);
#define Geom_temp_Defun(x) core::af_def(GeomPkg, #x, &geom__##x, ARGS_geom__##x, DECL_geom__##x, DOCS_geom__##x, __FILE__, __LINE__);
#define Units_temp_Defun(x) core::af_def(UnitsPkg, #x, &units__##x, ARGS_units__##x, DECL_units__##x, DOCS_units__##x, __FILE__, __LINE__);

#endif
