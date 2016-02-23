/*
    File: _AnchorRestraint_debugEvalSerialize.cc
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
// XML_EVAL
	node->attributeIfNotDefault("evalx1", this->eval.x1, 0.0 );
 	node->attributeIfNotDefault("evaly1", this->eval.y1, 0.0 );
 	node->attributeIfNotDefault("evalz1", this->eval.z1, 0.0 );
 	node->attributeIfNotDefault("evalxa", this->eval.xa, 0.0 );
 	node->attributeIfNotDefault("evalya", this->eval.ya, 0.0 );
 	node->attributeIfNotDefault("evalza", this->eval.za, 0.0 );
 	node->attributeIfNotDefault("evalka", this->eval.ka, 0.0 );
 	node->attributeIfNotDefault("evalEnergy", this->eval.Energy, 0.0 );
 	node->attributeIfNotDefault("evalfx1", this->eval.fx1, 0.0 );
 	node->attributeIfNotDefault("evalfy1", this->eval.fy1, 0.0 );
 	node->attributeIfNotDefault("evalfz1", this->eval.fz1, 0.0 );
 	node->attributeIfNotDefault("evaldhx1x1", this->eval.dhx1x1, 0.0 );
 	node->attributeIfNotDefault("evaldhy1y1", this->eval.dhy1y1, 0.0 );
 	node->attributeIfNotDefault("evaldhz1z1", this->eval.dhz1z1, 0.0 );
 	node->attributeIfNotDefault("evalohx1y1", this->eval.ohx1y1, 0.0 );
 	node->attributeIfNotDefault("evalohx1z1", this->eval.ohx1z1, 0.0 );
 	node->attributeIfNotDefault("evalohy1z1", this->eval.ohy1z1, 0.0 );
 	node->attributeIfNotDefault("evalAnchorDeviation", this->eval.AnchorDeviation, 0.0 );
