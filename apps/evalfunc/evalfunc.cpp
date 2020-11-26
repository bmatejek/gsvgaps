////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include "evalfunc.h"
#include "helper.cpp"
#include <limits>
#include <algorithm>

using namespace std;

// useful struct
struct Correlation {
  double score;
  int point;
};

////////////////////////////////////////////////////////////////////////
// Program Arguments
////////////////////////////////////////////////////////////////////////

// global macros
#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

// global program arguments
static int print_verbose;
static char *input_scene_name;
static char *object_name;

// global variables
static GSVScene *scene;
static int run_number = 0;
static vector<R3Point> truth_points = vector<R3Point>();
static vector<struct Correlation> correlations = vector<struct Correlation>();

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

// read in the GSVScene
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

// read in voc files
static int ReadVOCFiles(void) {
  // statistics
  RNTime start_time;
  start_time.Read();

  // get all truth points
  char truthFilename[128];
  snprintf(truthFilename, 127, "truth_data/%s/%02d/truth.txt", object_name, run_number);
  
  // open up the truth filename
  ifstream truth(truthFilename);
  if (truth.is_open()) {
    string line;
    while (getline(truth, line)) {
      vector<string> truthPoint = splitString(line, ',');
      
      double xcoord = atof(truthPoint[0].c_str());
      double ycoord = atof(truthPoint[1].c_str());
      double zcoord = atof(truthPoint[2].c_str());
      
      // add truth point
      truth_points.push_back(R3Point(xcoord, ycoord, zcoord));

    }
  } 
  else {
    fprintf(stderr, "Failed to open %s\n", truthFilename);
    return 0;
  }
  // close the file
  truth.close();

  // go through each run as input by user
  int gsv_run = ConvertIndividualRun(scene, run_number);
  GSVRun *run = scene->Run(gsv_run);
  // go through each segment in run
  for (int is = 0; is < run->NSegments(); is++) {
    GSVSegment *segment = run->Segment(is);
    // go through each panorama in segment
    for (int ip = 0; ip < segment->NPanoramas(); ip++) {
      GSVPanorama *panorama = segment->Panorama(ip);
      // go through each image in panorama
      for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
	GSVImage *image = panorama->Image(ii);
			
	char vocFilename[128];
	snprintf(vocFilename, 127, "voc_predictions/%s/%02d/%02d_%06d_%02d_UndistortedImage.txt", object_name, run_number, is, ip, ii);
	
	// open up the individual voc filename
	ifstream boxes(vocFilename);
	if (boxes.is_open()) {
	  string line;
	  // get parameters for each box in the file
	  while (getline(boxes, line)) {
	    vector<string> vocBox = splitString(line, ',');
	    
	    // ymin and ymax are switched (from 3 to 1 and vice versa) to take into consideration the changing coordinate system
	    double xmin = VOCToGSV(atof(vocBox[0].c_str()), true, GSV_HEIGHT, GSV_WIDTH);
	    double ymin = VOCToGSV(atof(vocBox[3].c_str()), false, GSV_HEIGHT, GSV_WIDTH);
	    double xmax = VOCToGSV(atof(vocBox[2].c_str()), true, GSV_HEIGHT, GSV_WIDTH);
	    double ymax = VOCToGSV(atof(vocBox[1].c_str()), false, GSV_HEIGHT, GSV_WIDTH);
	    R2Box vocPredictedBox = R2Box(xmin, ymin, xmax, ymax);
	    double score = atof(vocBox[5].c_str());
	    boolean found = false;
	    // go through every truth point to see if it fits in image
	    for (unsigned int it = 0; it < truth_points.size(); it++) {
	      R2Point imageCoord = image->UndistortedPosition(truth_points[it]);
	      if (imageCoord.X() > 0 && imageCoord.X() < GSV_WIDTH && imageCoord.Y() > 0 && imageCoord.Y() < GSV_HEIGHT) {
		
		// see if this image coordinate lies in the box range
		if (xmin < imageCoord.X() && imageCoord.X() < xmax && ymin < imageCoord.Y() && imageCoord.Y() < ymax) {
		  struct Correlation correlation;
		  correlation.score = score;
		  correlation.point = it;
		  
		  correlations.push_back(correlation);
		  found = true;
		  break;
		}
	      }
	    }
	    
	    // add to the list of boxes that have no correlation
	    if (!found) {
	      struct Correlation correlation;
	      correlation.score = score;
	      correlation.point = -1;
	      correlations.push_back(correlation);
	    }
	    
	  }
	}
	else {
	  fprintf(stderr, "Failed to open %s\n", vocFilename);
	  return 0;
	}
	// close the file
	boxes.close();

	// print the filename
	printf("Read %s...\n", vocFilename);
      }
    }
  }

  // print statistics
  if (print_verbose) {
    printf("Read VOC features\n");
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }

  // return OK status
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Determine Recall
////////////////////////////////////////////////////////////////////////

int determineRecall(void) {
  double bestScores[truth_points.size()];
  for (unsigned int i = 0; i < truth_points.size(); i++) {
    bestScores[i] = -1 * numeric_limits<double>::max();
  }
  // find the best scores for each point
  for (unsigned int i = 0; i < correlations.size(); i++) {
    double score = correlations[i].score;
    int point = correlations[i].point;
    if (point != -1 && score > bestScores[point]) {
      bestScores[point] = score;
    }
  }
  sort(bestScores, bestScores + truth_points.size());

  int wrongScores[truth_points.size()];
  for (unsigned int i = 0; i < truth_points.size(); i++) {
    wrongScores[i] = 0;
  }

  // go through all points and find false positives
  for (unsigned int i = 0; i < correlations.size(); i++) {
    double score = correlations[i].score;
    int point = correlations[i].point;

    if (point == -1) {
      for (unsigned int j = 0; j < truth_points.size(); j++) {
	if (score > bestScores[j]) {
	  wrongScores[j]++;
	}
	else {
	  break;
	}
      }
    }
  }

  char recallFilename[128];
  snprintf(recallFilename, 127, "evaluation_metrics/%s/%02d/voc_precision.txt", object_name, run_number);
  
  ofstream recallOutput;
  recallOutput.open(recallFilename);
  if (recallOutput.is_open()) {
    for (unsigned int i = 0; i < truth_points.size(); i++) {
      double recall = (i + 1) / (double) (i + wrongScores[i] + 1);
      recallOutput << recall;
      recallOutput << "\n";
    }
  }
  else {
    fprintf(stderr, "Failed to write to %s\n", recallFilename);
    return 0;
  }
  if (print_verbose)
    printf("Wrote %s\n", recallFilename);
  // return OK status
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Function
////////////////////////////////////////////////////////////////////////

void Usage(void) {
  fprintf(stderr, "evalfunc.exe [-v] [-input_scene] <input_scene> [-object] <object> -run_number <run_number>\n");
}

int parseArgs(int argc, char **argv) {
  // parse arguments
  argv++; argc--;
  
  while (argc > 0) {
    if (*argv[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-input_scene")) { argv++; argc--; input_scene_name = *argv; }
      else if (!strcmp(*argv, "-object")) { argv++; argc--; object_name = *argv; }
      else if (!strcmp(*argv, "-run_number")) { argv++; argc--; run_number = atoi(*argv); }
      else { fprintf(stderr, "Unrecognized argument: %s\n", *argv); return 0; }
    }
    else {
      if (!input_scene_name) { input_scene_name = *argv; }
      else if (!object_name) { object_name = *argv; }
    }
    argv++; argc--;
  }

  // return failure if there is no input scene
  if (!input_scene_name) {
    fprintf(stderr, "Need to include an input scene file\n");
    return 0;
  }

  // return failure if there is no object
  if (!object_name) {
    fprintf(stderr, "Need to include an object name\n");
    return 0;
  }

  // return OK status
  return 1;

}

////////////////////////////////////////////////////////////////////////
// Main 
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  
  // parse the command line arguments
  if (!parseArgs(argc, argv)) { Usage(); exit(-1); }

  // read in the scene
  scene = ReadScene(input_scene_name);
  if (!scene) { exit(-1); }
  
  // get precision and recall data
  if (!ReadVOCFiles()) { exit(-1); }

  // determine recall versus precision
  if (!determineRecall()) { exit(-1); }

  // return OK status 
  return 1;
}
