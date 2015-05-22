       
       
#ifndef	RandomGenerators_H
#define RandomGenerators_H

#include <clasp/core/common.h>
#include <boost/random.hpp>

#include <cando/chem/chemPackage.h>


namespace chem {

extern	boost::mt19937	_OnlyGenerator;

class Dumb_GaussianRandom;
    typedef	gctools::smart_ptr<Dumb_GaussianRandom>	RPGaussianRandom;
//class Dumb_GaussianRandom : public boost::enable_shared_from_this<Dumb_GaussianRandom>
class Dumb_GaussianRandom
{
private:
	boost::mt19937			*_Generator;
	boost::normal_distribution<>	*_Distribution;
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > *_VariateGenerator;

public:
	Dumb_GaussianRandom( double mean, double stdev );
	virtual ~Dumb_GaussianRandom();

	double	next();

};



// Constructor function

inline	RPGaussianRandom new_RPGaussianRandom(double mean, double stdev ) {
    return RPGaussianRandom(new Dumb_GaussianRandom(mean,stdev));
}

};
#endif
