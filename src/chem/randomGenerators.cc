/*
    File: randomGenerators.cc
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
#define	TURN_DEBUG_OFF


#include <cando/chem/randomGenerators.h>



namespace chem {
boost::mt19937	_OnlyGenerator;

};


#if 0
//
// Constructor
//
Dumb_GaussianRandom::Dumb_GaussianRandom( double mean, double stdev )
{
//    this->_Generator = new boost::mt19937();
    this->_Distribution = new boost::normal_distribution<> (mean,stdev);
    this->_VariateGenerator = new boost::variate_generator<boost::mt19937&,
	boost::normal_distribution<> > (_OnlyGenerator,*(this->_Distribution));
//    this->_VariateGenerator = new boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > (*(this->_Generator),*(this->_Distribution));
}

//
// Destructor
//
Dumb_GaussianRandom::~Dumb_GaussianRandom()
{
    delete this->_VariateGenerator;
    delete this->_Distribution;
//    delete this->_Generator;
}


double	Dumb_GaussianRandom::next()
{
    return (*(this->_VariateGenerator))();
}


#endif
