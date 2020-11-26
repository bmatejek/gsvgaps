#ifndef VOCFEATURE_H
#define VOCFEATURE_H

////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "GSV/GSV.h"
#include "R2Shapes/R2Shapes.h"
#include "R3Shapes/R3Shapes.h"
#include "RNBasics/RNBasics.h"
#include "R3CatmullRomSpline.h"
#include "poseoptimization.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

class VOCFeature {
 public:
  // constructors
  VOCFeature(int, int, int, int, double, double, double, double, double, GSVImage *);
  
  // accessors
  const R2Box BBox(void) const;
  const double Threshold(void) const;
  const GSVImage *Image(void) const;
  const double XMin(void) const;
  const double YMin(void) const;
  const double XMax(void) const;
  const double YMax(void) const;
  const int Run(void) const;
  const int Segment(void) const;
  const int Tapestry(void) const;
  const int Panorama(void) const;

 private:
  // instance variables
  R2Box bbox;
  double threshold;
  GSVImage *image;
  int run;
  int segment;
  int tapestry;
  int panorama;
};

class VOCFeatureList {
 public:
  // constructors
  VOCFeatureList();	
  
  // accessors/manipulators
  const std::vector<VOCFeature *> VOCFeatures(void) const;
  const VOCFeature *Feature(int) const;
  int ReadFile(const char *, GSVPoseOptimization *, double);
  int NFeatures();	
  
 private:
  void parseVOCFeatureLine(std::string, GSVPoseOptimization *, double);
  void addVOCFeature(VOCFeature *vocFeature);

  // instance variables
  std::vector<VOCFeature *> vocFeatures;
};

std::vector<std::string> splitString(std::string line, char delim);

#endif // VOCFEATURE_H