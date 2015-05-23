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
