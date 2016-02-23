/*
    File: otherPackageClasses.h
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
//
//

LISP_OTHER_PACKAGE(CorePkg,geom::HighlightedObject_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::Iterator_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::ConditionBase_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::ConditionWithStackTrace_O,core::ConditionBase_O)
LISP_OTHER_PACKAGE(CorePkg,core::CandoException_O,core::ConditionWithStackTrace_O)
LISP_OTHER_PACKAGE(CorePkg,core::LispError_O,core::CandoException_O)
LISP_OTHER_PACKAGE(CorePkg,geom::GraphicsObject_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,geom::Render_O,geom::GraphicsObject_O)
LISP_OTHER_PACKAGE(CorePkg,core::SequenceBaseVirtual_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::SequenceBase_O,core::SequenceBaseVirtual_O)
LISP_OTHER_PACKAGE(CorePkg,core::Sequence_O,core::SequenceBase_O)
