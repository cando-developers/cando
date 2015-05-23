#ifndef chem_units_H
#define chem_units_H

#if 0
#include <boost/mpl/list.hpp>
#include <boost/units/conversion.hpp>
#include <boost/units/io.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/prefixes.hpp>
using namespace boost::units;

#endif


namespace chem
{



    static const double Const_kJ_per_kCal = 4.1868;



inline double angstrom_to_nanometer(double val)
{
    return val/10.0;
}


inline double kCalPerAngstromSquared_to_kJPerNanometerSquared(double val)
{
	return val*Const_kJ_per_kCal*100.0;
}

inline double kJPerNanometerSquared_to_kCalPerAngstromSquared(double val)
{
    return val/100.0/Const_kJ_per_kCal;
}

inline double kCalPerRadianSquared_to_kJPerRadianSquared(double val)
{
    return val*Const_kJ_per_kCal;
}

inline double kJPerRadianSquared_to_kCalPerRadianSquared(double val)
{
    return val/Const_kJ_per_kCal;
}

inline double kCal_to_kJ(double val)
{
    return val*Const_kJ_per_kCal;
}

inline double kJ_to_kCal(double val)
{
    return val/Const_kJ_per_kCal;
}



};
#endif
