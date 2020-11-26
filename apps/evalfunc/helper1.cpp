////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

using namespace std;

#include "helper.h"

#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

////////////////////////////////////////////////////////////////////////
// Point2D Constructor
////////////////////////////////////////////////////////////////////////

Point2D ::
Point2D(R2Point point, coordinate_type system) :
	point(point),
	system(system)
{
}

////////////////////////////////////////////////////////////////////////
// Point2D Accessor Functions
////////////////////////////////////////////////////////////////////////

const R2Point Point2D ::
Point(void) const {
	return point;
}

const double Point2D ::
X(void) const {
	return point.X();
}

const double Point2D ::
Y(void) const {
	return point.Y();
}

const coordinate_type Point2D ::
System(void) const {
	return system;
}

////////////////////////////////////////////////////////////////////////
// Point2D Conversion Functions
////////////////////////////////////////////////////////////////////////

// convert to gsv coordinate system
void Point2D ::
ConvertToGSV(void) {
	if (system == PYGAME_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == VOC_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == IMG_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, GSV_WIDTH / 2 + point.X());
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 + point.Y());
	}
	system = GSV_COORDINATE_SYSTEM;
}

// convert to standard coordinate system
void Point2D ::
ConvertToStandard(void) {
	if (system == PYGAME_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == VOC_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == IMG_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, GSV_WIDTH / 2 + point.X());
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 + point.Y());
	}
	system = STANDARD_COORDINATE_SYSTEM;
}

// convert to pygame coordinate system
void Point2D ::
ConvertToPygame(void) {
	if (system == GSV_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == STANDARD_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == IMG_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, GSV_WIDTH / 2 + point.X());
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 - point.Y());
	}
	system = PYGAME_COORDINATE_SYSTEM;
}

// convert to voc coordinate system
void Point2D ::
ConvertToVOC(void) {
	if (system == GSV_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == STANDARD_COORDINATE_SYSTEM) {
		point.SetCoord(RN_Y, GSV_HEIGHT - point.Y());
	}
	else if (system == IMG_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, GSV_WIDTH / 2 + point.X());
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 - point.Y());
	}
	system = VOC_COORDINATE_SYSTEM;
}

// convert to image coordinate system
void Point2D ::
ConvertToIMG(void) {
	if (system == GSV_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, point.X() - GSV_WIDTH / 2);
		point.SetCoord(RN_Y, point.Y() - GSV_HEIGHT / 2);
	}
	else if (system == STANDARD_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, point.X() - GSV_WIDTH / 2);
		point.SetCoord(RN_Y, point.Y() - GSV_HEIGHT / 2);
	}
	else if (system == PYGAME_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, point.X() - GSV_WIDTH / 2);
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 - point.Y());
	}
	else if (system == VOC_COORDINATE_SYSTEM) {
		point.SetCoord(RN_X, point.X() - GSV_WIDTH / 2);
		point.SetCoord(RN_Y, GSV_HEIGHT / 2 - point.Y());
	}
	system = IMG_COORDINATE_SYSTEM;
}

////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////
