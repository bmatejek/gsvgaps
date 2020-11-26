#ifndef HELPER_H
#define HELPER_H

////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include "gaps.h"
#include "GSV/GSV.h"

////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

std::vector<std::string> splitString(std::string line, char delim);
int ConvertIndividualRun(GSVScene *scene, int run);
int IndividualRunToHumanReadable(GSVScene *scene, int run);
double VOCToImage(double coordinate, bool xCoord, int imageHeight, int imageWidth);
double VOCToGSV(double coordinate, bool xCoord, int imageHeight, int imageWidth);
double GSVToImage(double coordinate, bool xCoord, int imageHeight, int imageWidth);
double EstimatedImagePlaneDistance(double ymin, double ymax, double estimatedHeight, GSVImage *image);
double EstimatedImagePointDistance(double x, double y, double image_plane_distance, GSVImage *image);

#endif // HELPER_H
