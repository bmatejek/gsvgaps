////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <gaps.h>
#include <GSV/GSV.h>
#include "helper.h"

using namespace std;


////////////////////////////////////////////////////////////////////////
// Helper struct
////////////////////////////////////////////////////////////////////////

class VOCBox {
public:
  VOCBox(R2Box box, double threshold);

  R2Box box;
  double threshold;
};

VOCBox::
VOCBox(R2Box box, double threshold) :
  box(box),
  threshold(threshold)
{
}

typedef R2Box TruthBox;

class Image {
public:
  Image(GSVImage *Image);

  GSVImage *GSVImageReturn(void);
  vector<VOCBox *> VOCBoxes(void);
  vector<TruthBox *> TruthBoxes(void);

  void addVOCBox(VOCBox *);
  void addTruthBox(TruthBox *); 

  GSVImage *gsvImage;
  vector<VOCBox *> vocBoxes;
  vector<TruthBox *> truthBoxes;
};

Image::
Image(GSVImage *Image) :
  gsvImage(Image),
  vocBoxes(vector<VOCBox *>()),
  truthBoxes(vector<TruthBox *>())
{
}

GSVImage *Image::
GSVImageReturn(void) {
  return gsvImage;
}

vector<VOCBox *> Image::
VOCBoxes(void) {
  return vocBoxes;
}


vector <TruthBox *> Image::
TruthBoxes(void) {
  return truthBoxes;
}

void Image::
addVOCBox(VOCBox *box) {
  vocBoxes.push_back(box); 
}

void Image::
addTruthBox(TruthBox *box) {
  truthBoxes.push_back(box);
}

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static int run_number = 0;
static char *voc_directory = NULL;
static char *input_scene = NULL;
static GSVScene *scene = NULL;
static vector<Image *> images = vector<Image *>();

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

// read in the GSVScene from filename
static GSVScene *ReadScene(const char *filename) {
  // start statistics
  RNTime start_time;
  start_time.Read();
  
  // allocate google scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }
  
  // read scene
  if (!scene->ReadFile(filename, 0)) {
    delete scene;
    return NULL;
  }
  
  // print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
  // return scene
  return scene;
}

static int ReadVOCBoxes(void) {
  // HARD CODED IN  
  GSVRun *run = scene->Run(13);
  int run_number = 3;
  for (int is = 0; is < run->NSegments(); is++) {
    GSVSegment *segment = run->Segment(is);
    for (int ip = 0; ip < segment->NPanoramas(); ip++) {
      GSVPanorama *panorama = segment->Panorama(ip);
      for (int ii = 0; ii < panorama->NImages(); ii++) {
	// create new image
	GSVImage *gsvImage = panorama->Image(ii);
	Image *image = new Image(gsvImage);
	images.push_back(image);

	char vocFilename[128];
	snprintf(vocFilename, 127, "%s/%02d/%02d_%06d_%02d_UndistortedImage.txt", voc_directory, run_number, is, ip, ii);
	ifstream voc(vocFilename);
	if (voc.is_open()) {
	  string line;
	  while (getline(voc, line)) {
	    vector<string> args = splitString(line, ',');
	    
	    // parse box parameters
	    double xmin = atof(args[0].c_str());
	    double ymin = atof(args[1].c_str());
	    double xmax = atof(args[2].c_str());
	    double ymax = atof(args[3].c_str());
	    double score = atof(args[5].c_str());

	    R2Box box = R2Box(xmin, ymin, xmax, ymax);
	    VOCBox *vocBox = new VOCBox(box, score);
	    image->addVOCBox(vocBox);
	  }
	  

	}
	else {
	  fprintf(stderr, "Error reading file %s\n", vocFilename);
	}
      }
    }
  }

  // RETURN SUCCESS
  return 1;
}   


////////////////////////////////////////////////////////////////////////
// Comparing Functions
////////////////////////////////////////////////////////////////////////

static double IntersectionOverUnion(R2Box box1, R2Box box2) {
  /// create intersection box
  R2Box intersect = R2Box(box1);
  intersect.Intersect(box2);
  
  double areaIntersect = intersect.Area();
  double areaUnion = box1.Area() + box2.Area() - intersect.Area();

  return areaIntersect / areaUnion;
}


////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static void Usage(void) {
  /* TODO THIS */
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
  // parse args
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene = *argv; }
      else if (!strcmp(*argv, "-voc_directory")) { argc--; argv++; voc_directory = *argv; }
      else if (!strcmp(*argv, "-run_number")) { argc--; argv++; run_number = atoi(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!input_scene) { input_scene = *argv; }
      else if (!voc_directory) { voc_directory = *argv; }
    }
    argc--; argv++;
  }
  
  if (!voc_directory) { fprintf(stderr, "No voc directory given.\n"); return 0; }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0; }

  scene = ReadScene(input_scene);
  if (!scene) { return 0; }

  if (!ReadVOCBoxes()) { return 0;}

  // return SUCCESS
  return 1;
}
