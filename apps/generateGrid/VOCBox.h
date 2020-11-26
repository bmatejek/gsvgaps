#ifndef VOC_BOX_H
#define VOC_BOX_H

////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"
#include "generateGrid.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

class VOCBox {
 public:
  // constructor
  VOCBox(double xmin, double ymin, double xmax, double ymax, double score, GSVImage *Image, ObjectOfInterest *object);
  
  // accessors/manipulators
  const double XMin(void) const;
  const double YMin(void) const;
  const double XMax(void) const;
  const double YMax(void) const;
  const double Score(void) const;
  const GSVImage *Image(void) const;
  const ObjectOfInterest *Object(void) const;
  double Width(void);
  double Height(void);
  R2Point Center(void);
  double CenterX(void);
  double CenterY(void);
  double ImagePlaneDistance(void);
  R3Point GlobalPoint(void);

 private:
  // instance variables
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  double score;
  GSVImage *image;
  ObjectOfInterest *object;
};

#endif // VOC_BOX_H
