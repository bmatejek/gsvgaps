////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

// transform a voc coordinate into an image coordinate
double VOCToImage(double coordinate, bool xCoord, int imageHeight, int imageWidth) {
  if (xCoord) {
    return coordinate - imageWidth / 2;
  }
  else {
    return imageHeight / 2 - coordinate;
  }
}

// transform a voc coordinate into a gsv coordinate
double VOCToGSV(double coordinate, bool xCoord, int imageHeight, int imageWidth) {
  if (xCoord) {
    return coordinate;
  }
  else {
    return imageHeight - coordinate;
  }
}

// split the line by delim and return the vector
vector<string> splitString(string line, char delim) {
  vector<string> arguments;
  int length = line.length();
  int string_start = 0;
  
  // split the string by delim
  for (int i = 0; i < length; i++) {
    if (line.at(i) == delim) {
      int string_end = i;
      if (string_start != string_end) {
	arguments.push_back(line.substr(string_start, string_end - string_start));
      }
      string_start = i + 1;
    }
  }
  
  // add the last substring
  if (string_start != length)
    arguments.push_back(line.substr(string_start));
  
  // return the number of arguments found
  return arguments;
}
