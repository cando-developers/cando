/* @(#)dimension.fwd.h
 */

#ifndef _units_Dimension_fwd_H
#define _units_Dimension_fwd_H


namespace units
{

	typedef enum {distanceMeter=0,
		      weightKilogram=1,
		      timeSecond=2,
		      currentAmpere=3,
		      temperatureKelvin=4,
		      luminosityCandela=5,
		      amountMole=6,
		      radianAngle=7,
		      NumBaseDimensions=8
	} BaseDimensionEnum;





FORWARD(Dimension);



}; /* units */

#endif /* _units_Dimension_fwd_H */
