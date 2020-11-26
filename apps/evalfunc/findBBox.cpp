////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

using namespace std;

#include "findBBox.h"

////////////////////////////////////////////////////////////////////////
// Program Arguments
////////////////////////////////////////////////////////////////////////

// global macros
#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

int findBBox(CorrelateImage *correlation, int primaryImage, int secondaryImage) {
	R3Point globalObject = PointFromImage(correlation, primaryImage, 163, R2Point(482.0, 1245.0));
	
	Image imageData = correlation->Images(secondaryImage);
	GSVImage *image = imageData.image;
	
	// get camera parameters
	GSVCamera *camera = image->Camera();
	RNAngle xfov = 0.5 * camera->XFov();
	RNAngle yfov = 0.5 * camera->YFov();
	GSVPose pose = image->Pose();
	R3Point viewpoint = pose.Viewpoint();
	R3Vector towards = pose.Towards();
	R3Vector up = pose.Up();
	R3Vector right = pose.Right();
	
	// get vector from viewpoint to object
	R3Vector viewpointToObject = globalObject - viewpoint;
	double image_plane_distance = viewpointToObject.Dot(towards);
	double offsetX = viewpointToObject.Dot(right);
	double offsetY = viewpointToObject.Dot(up);
	printf("Vector offset: %f %f\n", offsetX, offsetY);
	// find out how long the vectors to the edge of the image are
	R3Vector dx = right * image_plane_distance * tan(xfov);
	R3Vector dy = up * image_plane_distance * tan(yfov);
	double xDist = dx.Length();
	double yDist = dy.Length();
	printf("Distance (dx, dy): %f %f\n", xDist, yDist);
	// get pixel offset
	double pixelX = offsetX / xDist * GSV_WIDTH / 2;
	double pixelY = offsetY / yDist * GSV_HEIGHT / 2;
	printf("Pixel Offset: %f %f\n", pixelX, pixelY);
	Point2D imgPoint = Point2D(R2Point(pixelX, pixelY), IMG_COORDINATE_SYSTEM);
	imgPoint.ConvertToVOC();
	
	// get object height 
	double height  = (3.0 / 8.0) * GSV_HEIGHT / (image_plane_distance * tan(yfov));
	printf("%f\n", height);
	printf("VOC Image Offset: %f %f\n", imgPoint.X(), imgPoint.Y());
	
	printf("Image Plane Distance: %f\n", image_plane_distance);
	
	// return OK status
	return 1;
}

R3Point PointFromImage(CorrelateImage *correlation, int imageIndex, double height, R2Point center) {
	// find global variables
	Image imageData = correlation->Images(imageIndex);
	GSVImage *image = imageData.image;
	VOCFeature *vocFeature = imageData.vocFeatures[0];
	
	printf("%d %d %d %d\n", image->RunIndex(), image->SegmentIndex(), image->TapestryIndex(), image->PanoramaIndex());
	
	// get camera parameters
	GSVCamera *camera = image->Camera();
	RNAngle xfov = 0.5 * camera->XFov();
	RNAngle yfov = 0.5 * camera->YFov();
	GSVPose pose = image->Pose();
	R3Point viewpoint = pose.Viewpoint();
	R3Vector towards = pose.Towards();
	R3Vector up = pose.Up();
	R3Vector right = pose.Right();
	
	// find  the image plane distance if the picture were on the stop sign
	//double height = height; //163.0; // in pixels
	double image_plane_distance = (3.0 / 8.0) * ((double) GSV_HEIGHT) / (height * tan(yfov));
	
	// relVOC is the center of the bounding box - relIMG is the same point relative to image coordinates
	// find coordinates relative to central image origin
	Point2D point = Point2D(center, VOC_COORDINATE_SYSTEM); //Point2D(R2Point(482.0, 1245.0), VOC_COORDINATE_SYSTEM);
	point.ConvertToIMG();
	
	// get vector to the ends of the image
	R3Vector dx = right * image_plane_distance * tan(xfov);
	R3Vector dy = up * image_plane_distance * tan(yfov);
	
	// get vectors from image origin by scaling dx, dy
	double maxX = GSV_WIDTH / 2;
	double maxY = GSV_HEIGHT / 2;
	dx = dx * point.X() / maxX;
	dy = dy * point.Y() / maxY;
	R3Vector fromImageOrigin = dx + dy;
	R3Vector fromCamera = towards * image_plane_distance + fromImageOrigin;
	R3Point objectGlobal = viewpoint + fromCamera;
	return objectGlobal;
}

double MatchingBBoxes(CorrelateImage *correlation, int primaryIndex, int secondaryIndex) {
	R3Point primaryPoint = PointFromImage(correlation, primaryIndex, 163, R2Point(482.0, 1245.0));
	R3Point secondaryPoint = PointFromImage(correlation, secondaryIndex, 139, R2Point(208, 1250.0));
	printf("First point: %f %f %f\n", primaryPoint.X(), primaryPoint.Y(), primaryPoint.Z());
	printf("Second point: %f %f %f\n", secondaryPoint.X(), secondaryPoint.Y(), secondaryPoint.Z());
	double distance = R3Distance(primaryPoint, secondaryPoint);
	// create some probability based on distance, bbox error, etc.
	return 0.5;
}