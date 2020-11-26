////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include "vocFeature.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// VOCFeature Constructors
////////////////////////////////////////////////////////////////////////

VOCFeature ::
VOCFeature(int run, int segment, int tapestry, int panorama, double xmin, double ymin, double xmax, double ymax, double threshold, GSVImage *image) :
  bbox(R2Box(xmin, ymin, xmax, ymax)),
  threshold(threshold),
  image(image),
  run(run),
  segment(segment),
  tapestry(tapestry),
  panorama(panorama)
{
}

////////////////////////////////////////////////////////////////////////
// VOCFeature Accessors
////////////////////////////////////////////////////////////////////////

const R2Box VOCFeature ::
BBox(void) const {
  return bbox;
}

const double VOCFeature ::
Threshold(void) const {
  return threshold;
}

const GSVImage *VOCFeature ::
Image(void) const {
  return image;
}

const double VOCFeature ::
XMin(void) const {
  return bbox.XMin();
}

const double VOCFeature ::
YMin(void) const {
  return bbox.YMin();
}

const double VOCFeature ::
XMax(void) const {
  return bbox.XMax();
}

const double VOCFeature ::
YMax(void) const {
  return bbox.YMax();
}

const int VOCFeature ::
Run(void) const {
  return run;
}

const int VOCFeature ::
Segment(void) const {
  return segment;
}

const int VOCFeature ::
Tapestry(void) const {
  return tapestry;
}

const int VOCFeature ::
Panorama(void) const {
  return panorama;
}

////////////////////////////////////////////////////////////////////////
// VOCFeatureList Constructors
////////////////////////////////////////////////////////////////////////

VOCFeatureList ::
VOCFeatureList(void) {
  vocFeatures = vector<VOCFeature *>();
}

////////////////////////////////////////////////////////////////////////
// VOCFeatureList Accessors/Manipulators
////////////////////////////////////////////////////////////////////////

const vector<VOCFeature *> VOCFeatureList ::
VOCFeatures(void) const {
  return vocFeatures;
}

void VOCFeatureList ::
addVOCFeature(VOCFeature *vocFeature) {
  vocFeatures.push_back(vocFeature);
}

int VOCFeatureList ::
NFeatures() {
  return vocFeatures.size();
}

const VOCFeature *VOCFeatureList ::
Feature(int i) const {
  return vocFeatures[i];
}

////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////

// convert run from normal indices (1, 2, 3) to gsv indices (1, 10, 11 ... 19, 2, ...)
static int convertRun(int run) {
  int run_number = run;
  if (run >= 2 && run <= 9) run_number += 10;
  else if (run >= 10 && run <= 19) run_number -= 8;
  return run_number;
}

static int indexFromIndices(GSVScene *scene, int run_number, int segment_number, int tapestry_index, int panorama_index) {
	run_number = convertRun(run_number);
	int offset = 0;
	for (int i = 0; i < run_number; i++) {
		offset += scene->Run(i)->NImages();
	}
	for (int i = 0; i < segment_number; i++) {
		offset += scene->Run(run_number)->Segment(segment_number)->NImages();
	}
	return 9 * tapestry_index + panorama_index + offset;
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

////////////////////////////////////////////////////////////////////////
// Ascii Input and Output
////////////////////////////////////////////////////////////////////////

// given the line create a corresponding VOCImageFeature object
void VOCFeatureList ::
parseVOCFeatureLine(string line, GSVPoseOptimization *optimization, double threshold) {
  vector<string> arguments = splitString(line, '*');
	
  // get segment number, tapestry number, and panorama index for all	
  int run_number = atoi(arguments[0].c_str());
  int segment_number = atoi(arguments[1].c_str());
  int tapestry_index = atoi(arguments[2].c_str());
  int panorama_index = atoi(arguments[3].c_str());
  
  // get image
  int imageIndex = indexFromIndices(optimization->scene, run_number, segment_number, tapestry_index, panorama_index);	
  GSVImage *image = optimization->images[imageIndex];
  
  // find all the VOC features in this image
  for (unsigned int i = 4; i < arguments.size(); i += 4) {
    double x1 = atof(arguments[i].c_str());
    double y1 = atof(arguments[i + 1].c_str());
    double x2 = atof(arguments[i + 2].c_str());
    double y2 = atof(arguments[i + 3].c_str());
    
    // create VOCImageFeature
    VOCFeature *vocFeature = new VOCFeature(run_number, segment_number, tapestry_index, panorama_index, x1, y1, x2, y2, threshold, image);
    addVOCFeature(vocFeature);
  }
}

int VOCFeatureList ::
ReadFile(const char *filename, GSVPoseOptimization *optimization, double threshold) {
  // open file
  ifstream fp;
  fp.open(filename, ios::in);
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }
  
  // read through each line in the file
  for (string line; getline(fp, line); ) {
    parseVOCFeatureLine(line, optimization, threshold);
  }
  
  // close file
  fp.close();
  
  // return success
  return 1;
}