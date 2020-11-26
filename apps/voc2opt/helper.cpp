////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

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


int ConvertIndividualRun(GSVScene *scene, int run) {
  int nruns = scene->NRuns();
  // function does not support more than 100 runs
  if (nruns > 100) {
    fprintf(stderr, "%s on line %d of file %s assumes the scene has fewer than 100 runs (0 to 99).", __func__, __LINE__, __FILE__);
    exit(-1);
  }
  // convert the number
  if (run == 0) {
    return 0;
  }
  // case if run is less than 10
  else if (run < 10) {
    int numberBehind = 0;
    int mult10 = 10 * run;
    if (mult10 < nruns)
      numberBehind += nruns - mult10;
    numberBehind += 9 - run;
    return nruns - 1 - numberBehind;
  }
  // case if run is greater than 10
  else {
    int leadingDigit = run / 10;
    int trailingDigit = run % 10;
    return 11 * leadingDigit - 9 + trailingDigit;
  }
}

int IndividualRunToHumanReadable(GSVScene *scene, int run) {
  // convert run into human readable format
  if (run == 0)
    return 0;

  for (int i = 1; i < scene->NRuns(); i++) {
    if (ConvertIndividualRun(scene, i) == run) return i;
  }
  
  return -1;
}
