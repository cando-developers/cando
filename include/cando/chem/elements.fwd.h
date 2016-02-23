/*
    File: elements.fwd.h
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
#ifndef elements_FWD_H
#define elements_FWD_H


#include <cando/chem/chemPackage.h>


namespace chem
{




	/*! Enough room for all elements and then some special
	 * element entries like LonePairs, Deuterium, UndefinedElement, Dummy element
	 */

typedef enum { 
element_Dummy=0,
element_LP=200,
element_D=201,
element_T=202,
element_Undefined=203,
element_None=204,
element_Wildcard=205,
element_MAX=250,
element_Ac=89,
element_Ag=47,
element_Al=13,
element_Am=95,
element_Ar=18,
element_As=33,
element_At=85,
element_Au=79,
element_B=5,
element_Ba=56,
element_Be=4,
element_Bh=107,
element_Bi=83,
element_Bk=97,
element_Br=35,
element_C=6,
element_Ca=20,
element_Cd=48,
element_Ce=58,
element_Cf=98,
element_Cl=17,
element_Cm=96,
element_Co=27,
element_Cr=24,
element_Cs=55,
element_Cu=29,
element_Db=105,
element_Ds=110,
element_Dy=66,
element_Er=68,
element_Es=99,
element_Eu=63,
element_F=9,
element_Fe=26,
element_Fm=100,
element_Fr=87,
element_Ga=31,
element_Gd=64,
element_Ge=32,
element_H=1,
element_He=2,
element_Hf=72,
element_Hg=80,
element_Ho=67,
element_Hs=108,
element_I=53,
element_In=49,
element_Ir=77,
element_K=19,
element_Kr=36,
element_La=57,
element_Li=3,
element_Lr=103,
element_Lu=71,
element_Md=101,
element_Mg=12,
element_Mn=25,
element_Mo=42,
element_Mt=109,
element_N=7,
element_Na=11,
element_Nb=41,
element_Nd=60,
element_Ne=10,
element_Ni=28,
element_No=102,
element_Np=93,
element_O=8,
element_Os=76,
element_P=15,
element_Pa=91,
element_Pb=82,
element_Pd=46,
element_Pm=61,
element_Po=84,
element_Pr=59,
element_Pt=78,
element_Pu=94,
element_Ra=88,
element_Rb=37,
element_Re=75,
element_Rf=104,
element_Rg=111,
element_Rh=45,
element_Rn=86,
element_Ru=44,
element_S=16,
element_Sb=51,
element_Sc=21,
element_Se=34,
element_Sg=106,
element_Si=14,
element_Sm=62,
element_Sn=50,
element_Sr=38,
element_Ta=73,
element_Tb=65,
element_Tc=43,
element_Te=52,
element_Th=90,
element_Ti=22,
element_Tl=81,
element_Tm=69,
element_U=92,
element_V=23,
element_W=74,
element_Xe=54,
element_Y=39,
element_Yb=70,
element_Zn=30,
element_Zr=40 
} Element;



typedef	enum { hybridization_undefined=0,
	       hybridization_sp=1,
	       hybridization_sp2=2,
	       hybridization_sp3=3,
	       hybridization_MAX=4} Hybridization;



Element elementFromSymbol(core::Symbol_sp sym);





};

#endif // elements_FWD_H
