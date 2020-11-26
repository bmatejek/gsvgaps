////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"
#include <dirent.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static int blur_grid = 0;
static int scale_factor = 1;
static char *voc_filename = NULL;
static char *prob_filename = NULL;

enum {
  GSV_WIDTH = 1936,
  GSV_HEIGHT = 2592,
};

////////////////////////////////////////////////////////////////////////
// Create Probability Files
////////////////////////////////////////////////////////////////////////

static int CreateProbabilityGrid() {
  // create grid
  R2Grid *grid = new R2Grid(GSV_WIDTH, GSV_HEIGHT);
  if (!grid) { fprintf(stderr, "Failed to allocate memory for grid.\n"); return 0; }
      
  char fullPathName[128];
  snprintf(fullPathName, 128, "%s", voc_filename);
  ifstream vocFile;
  vocFile.open(fullPathName);
  if (vocFile.is_open()) {
    string line;
    while (getline(vocFile, line)) {
      vector<string> args = splitString(line, ',');
      
      if (args.size() != 6) {
	fprintf(stderr, "Failed to read %s\n", fullPathName);
	fprintf(stderr, "Line has too many values\n");
	return 0;
      }
      
      // get parameters for voc box
      double xmin = atof(args[0].c_str());
      double ymin = atof(args[1].c_str());
      double xmax = atof(args[2].c_str());
      double ymax = atof(args[3].c_str());
      
      // normalize the score between 0 and 1
      double score;
      // create the probability values by only considering scores better than 0.5
      if (score > 1.0) {
	score = 1.0;
      } else if (score > -0.5) {
	score = (atof(args[5].c_str()) + 0.5) / 1.5;
      } else {
	score = 0.0;
      }
      
      // add to the score of the grid
      for (int x = (int) floor(xmin); x < (int) ceil(xmax); x++) {
	for (int y = (int) floor(ymin); y < (int) ceil(ymax); y++) {
	  int newXCoord = (int) round(VOCToGSV((double) x, true, GSV_HEIGHT, GSV_WIDTH));
	  int newYCoord = (int) round(VOCToGSV((double) y, false, GSV_HEIGHT, GSV_WIDTH));
	  if (score > grid->GridValue(newXCoord, newYCoord)) {
	    grid->SetGridValue(newXCoord, newYCoord, score);
	  }
	}
      }
    }
    if (print_verbose)  {
      printf("Reading %s\n", fullPathName);
    }
  }
  else {
    fprintf(stderr, "Failed to open %s\n", fullPathName);
    // return FAILURE
    return 0;
  }
  
  // blur the grid if desired
  if (blur_grid) {
    // blur the grid by twice the scale_factor
    grid->Blur(2 * scale_factor);
  }
  
  // save the grid to the accompanying probability file 
  char probabilityPathName[128];
  snprintf(probabilityPathName, 128, "%s", prob_filename);

  R2Grid *scaleGrid = new R2Grid(GSV_WIDTH / scale_factor, GSV_HEIGHT / scale_factor);
  for (int x = 0; x < scaleGrid->XResolution(); x++) {
    for (int y = 0; y < scaleGrid->YResolution(); y++) {
      scaleGrid->SetGridValue(x, y, grid->GridValue(x * scale_factor, y * scale_factor));
    }
  }
  
  // write grid to grd file
  scaleGrid->WriteFile(probabilityPathName);
  fprintf(stderr, "Writing %s\n", probabilityPathName);
  
  // delete the grids
  delete grid;
  delete scaleGrid;
  
  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static void Usage(void) {
  fprintf(stderr, "voc2gsv <voc_filename> <probability_filename> -scale_factor <scale_factor> -blur_grid -v\n");
  fprintf(stderr, "To run multiple instances: \n");
  fprintf(stderr, "python voc2gsvScript.py voc_directory/ probability_directory/ <scale_factor> <?blur_grid>\n");
  fprintf(stderr, "chmod 700 voc2gsv_executable\n");
  fprintf(stderr, "./voc2gsv_executable\n");
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
  // parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-blur_grid")) { blur_grid = 1; }
      else if (!strcmp(*argv, "-scale_factor")) { argc--; argv++; scale_factor = atoi(*argv); }
      else if (!strcmp(*argv, "-usage")) { return 0; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!voc_filename) { voc_filename = *argv; }
      else if (!prob_filename) { prob_filename = *argv; }
    }
    argc--; argv++;
  }

  if (!voc_filename) { fprintf(stderr, "Need voc directory\n"); return 0; }
  if (!prob_filename) { fprintf(stderr, "Need probability directory\n"); return 0; }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0;}
  
  // Create Probability files
  if (!CreateProbabilityGrid()) { return 0; }

  // return SUCCESS
  return 1;
}
