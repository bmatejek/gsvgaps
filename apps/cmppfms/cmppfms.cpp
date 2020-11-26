////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <gaps.h>
#include <GSV/GSV.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static char *truth_grid = NULL;
static char *test_grid = NULL;

////////////////////////////////////////////////////////////////////////
// Grid Comparing Function
////////////////////////////////////////////////////////////////////////

struct ScorePair {
  int x_coordinate;
  int y_coordinate;
  double score;
};

static bool ScoreSorter(ScorePair const& a, ScorePair const& b) {
  return a.score > b.score;
}

static int CompareGrids(void) {
  // open up the two grids
  R2Grid *truth = new R2Grid();
  truth->ReadFile(truth_grid);
  if (!truth) {
    fprintf(stderr, "Failed to read %s\n", truth_grid);
    return 0;
  }

  R2Grid *test = new R2Grid();
  test->ReadFile(test_grid);
  if (!test) {
    fprintf(stderr, "Failed to read %s\n", test_grid);
    return 0;
  }

  if (test->NEntries() != truth->NEntries()) {
    fprintf(stderr, "Resolutions of grids do not match\n");
    return 0;
  }
  
  // create a vector of real objects
  vector<int> x_coordinates = vector<int>();
  vector<int> y_coordinates = vector<int>();
  for (int x = 0; x < truth->XResolution(); x++) {
    for (int y = 0; y < truth->YResolution(); y++) {
      if (truth->GridValue(x, y) > 0) {
	x_coordinates.push_back(x);
	y_coordinates.push_back(y);
      }
    }
  }
  
  // create array for test
  int resolution = test->NEntries();
  ScorePair *scores = (ScorePair *) malloc(sizeof(ScorePair) * resolution);
  int i = 0;
  for (int x = 0; x < test->XResolution(); x++) {
    for (int y = 0; y < test->YResolution(); y++, i++) {
      scores[i].x_coordinate = x;
      scores[i].y_coordinate = y;
      scores[i].score = test->GridValue(x, y);
    }
  }
  
  // sort the scores
  sort(scores, scores+resolution, &ScoreSorter);

  int matches_found = 0;
  int false_positives = 0;
  for (int i = 0; i < resolution; i++) {
    bool match = false;
    for (unsigned int j = 0; j < x_coordinates.size(); j++) {
      int distanceSquared = (x_coordinates[j] - scores[i].x_coordinate) * (x_coordinates[j] - scores[i].x_coordinate) + (y_coordinates[j] - scores[i].y_coordinate) * (y_coordinates[j] - scores[i].y_coordinate);
      if (distanceSquared <= 16) {
	match = true;
	break;
      }
    }
    if (match) {
      matches_found++;
      
      printf("%f\n", matches_found / ((double) (matches_found + false_positives)));
    }
    // if match not found precision decreases
    else {
      false_positives++;
    }
    // found all matches
    if (matches_found == (int) x_coordinates.size())
      break;
  }
  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static void Usage(void) {
  fprintf(stderr, "cmppfms [-v] <truth_grid> <compare_grid>\n");
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
  // parse args
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!truth_grid) { truth_grid = *argv; }
      else if (!test_grid) { test_grid = *argv; }
    }
    argc--; argv++;
  }
 
  if (!truth_grid) { fprintf(stderr, "Need to include truth grid\n"); return 0; }
  if (!test_grid) { fprintf(stderr, "Neet to include test grid\n"); return 0; }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0; }

  if (!CompareGrids()) { fprintf(stderr, "Failed to output precision and recall\n"); return 0; }

  // return SUCCESS
  return 1;
}
