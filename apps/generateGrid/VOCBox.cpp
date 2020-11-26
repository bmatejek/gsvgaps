////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "VOCBox.h"

////////////////////////////////////////////////////////////////////////
// VOCBox class functions
////////////////////////////////////////////////////////////////////////

// constructor
VOCBox ::
VOCBox(double xmin, double ymin, double xmax, double ymax, double score, GSVImage *image, ObjectOfInterest *object) :
  xmin(xmin),
  ymin(ymin),
  xmax(xmax),
  ymax(ymax),
  score(score),
  image(image),
  object(object)
{
  xmin = VOCToImage(xmin, true, image->Height(), image->Width());
  ymin = VOCToImage(ymin, false, image->Height(), image->Width());
  xmax = VOCToImage(xmax, true, image->Height(), image->Width());
  ymax = VOCToImage(ymax, false, image->Height(), image->Width());
}

// accessors/manipulators
const double VOCBox ::
XMin(void) const {
  return xmin;
}

const double VOCBox ::
YMin(void) const {
  return ymin;
}

const double VOCBox ::
XMax(void) const {
  return xmax;
}

const double VOCBox ::
YMax(void) const {
  return ymax;
}

const double VOCBox ::
Score(void) const {
  return score;
}

const GSVImage *VOCBox ::
Image(void) const {
  return image;
}

const ObjectOfInterest *VOCBox ::
Object(void) const {
  return object;
}

double VOCBox ::
Width(void) {
  return xmax - xmin;
}

double VOCBox ::
Height(void) {
  return ymax - ymin;
}

R2Point VOCBox ::
Center(void) {
  return R2Point(CenterX(), CenterY());
}

double VOCBox ::
CenterX(void) {
  return (xmax + xmin) / 2;
}

double VOCBox ::
CenterY(void) {
  return (ymax + ymin) / 2;
}

double VOCBox ::
ImagePlaneDistance(void) {
  // get camera parameters
  GSVCamera *camera = image->Camera();
  RNAngle yfov = 0.5 * camera->YFov();
  
  // get the object height
  double objectHeight = object->Height();
  // get image plane distance and return it
  double image_plane_distance = objectHeight / 2 * ((double) image->Height()) / (Height() * tan(yfov));
  return image_plane_distance;
}

R3Point VOCBox ::
GlobalPoint(void) {
  // get camera parameters
  GSVCamera *camera = image->Camera();
  RNAngle xfov = 0.5 * camera->XFov();
  RNAngle yfov = 0.5 * camera->YFov();
  GSVPose pose = image->Pose();
  R3Point viewpoint = pose.Viewpoint();
  R3Vector towards = pose.Towards();
  R3Vector up = pose.Up();
  R3Vector right = pose.Right();
  
  // get image plane
  double image_plane_distance = ImagePlaneDistance();
  
  // get vectors on image
  R3Vector dx = right * image_plane_distance * tan(xfov);
  R3Vector dy = up * image_plane_distance * tan(yfov);
  
  // get vectors from image origin by scaling dx, dy
  double maxX = image->Width() / 2;
  double maxY = image->Height() / 2;
  dx = dx * CenterX() / maxX;
  dy = dy * CenterY() / maxY;
  R3Vector fromImageOrigin = dx + dy;
  R3Vector fromCamera = towards * image_plane_distance + fromImageOrigin;
  R3Point objectGlobal = viewpoint + fromCamera;
  return objectGlobal;
}
