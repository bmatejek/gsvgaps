#ifndef HELPER_H
#define HELPER_H

////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include "gaps.h"
#include "GSV/GSV.h"
#include "R3CatmullRomSpline.h"
#include "poseoptimization.h"

////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

std::vector<std::string> splitString(std::string line, char delim);
int ConvertIndividualRun(GSVScene *scene, int run);
int IndividualRunToHumanReadable(GSVScene *scene, int run);
double VOCToImage(double coordinate, bool xCoord, int imageHeight, int imageWidth);

#endif // HELPER_H
