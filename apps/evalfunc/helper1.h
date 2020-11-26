#ifndef HELPER_H
#define HELPER_H

////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include "R2Shapes/R2Shapes.h"
#include "R3Shapes/R3Shapes.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

typedef enum {
	GSV_COORDINATE_SYSTEM,
	STANDARD_COORDINATE_SYSTEM,
	PYGAME_COORDINATE_SYSTEM,
	VOC_COORDINATE_SYSTEM,
	IMG_COORDINATE_SYSTEM,
} coordinate_type;

// struct for point in various coordinate systems
class Point2D {
public:
	Point2D(R2Point point, coordinate_type system);
	
	// accessor functions
	const R2Point Point(void) const;
	const double X(void) const;
	const double Y(void) const;
	const coordinate_type System(void) const;
	
	// conversion functions
	void ConvertToGSV(void);
	void ConvertToStandard(void);
	void ConvertToPygame(void);
	void ConvertToVOC(void);
	void ConvertToIMG(void);
	
private:
	// instance variables
	R2Point point;
	coordinate_type system; 
};


////////////////////////////////////////////////////////////////////////
// Visible Functions
////////////////////////////////////////////////////////////////////////

#endif // HELPER_H