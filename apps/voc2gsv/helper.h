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
double VOCToImage(double coordinate, bool xCoord, int imageHeight, int imageWidth);
double VOCToGSV(double coordinate, bool xCoord, int imageHeight, int imageWidth);

#endif // HELPER_H
