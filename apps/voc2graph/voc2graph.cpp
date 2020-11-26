////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <dirent.h>
#include "gaps.h"
#include "GSV/GSV.h"
#include "gco/gco.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static int run_number = 0;
static int scale_factor = 0;
static char *input_scene = NULL;
static char *voc_directory = NULL;
static char *input_grid = NULL;
static char *output_filename = NULL;

// program global variables
static int number_images = 0;
static GSVScene *scene = NULL;
static R3Grid *probability_grid = NULL;

#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

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

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

// read in the GSVScene from filename
static int ReadScene() {
  // start statistics
  RNTime start_time;
  start_time.Read();
  
  // allocate google scene
  scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return 0;
  }
 
  // read scene
  if (!scene->ReadFile(input_scene, 0)) {
    delete scene;
    return 0;
  }
  
  // print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", input_scene);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
  // return SUCCESS
  return 1;
}

static int ReadGrid() {
  RNTime start_time;
  start_time.Read();

  probability_grid = new R3Grid();
  if (!probability_grid) {
    fprintf(stderr, "Unable to allocate grid memory\n");
    return 0;
  }

  if (!probability_grid->ReadFile(input_grid)) {
    delete probability_grid;
    fprintf(stderr, "Unable to read grid from %s\n", input_grid);
    return 0;
  }

  
  if (print_verbose) {
    printf("Read grid %s ...\n", input_grid);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }

  return 1;
}

// write grid to filename
static void WriteGrid(const char *filename, R3Grid *grid) {
  RNTime start_time;
  start_time.Read();
  
  // write the grid
  if (!grid->WriteGridFile(filename)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return;
  }

  // print statistics
  if (print_verbose) {
    printf("Write grid to %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
}


////////////////////////////////////////////////////////////////////////
// Create GraphCut Input
////////////////////////////////////////////////////////////////////////

static int SiteIDNumberFromImage(int number_image, int x, int y) {
  int sites_per_image = (int) (ceil(GSV_WIDTH / (double) scale_factor) * ceil(GSV_HEIGHT / (double) scale_factor));
  int sites_before_other_images = sites_per_image * number_image;
  int sites_before_within_image = x + GSV_WIDTH * y;
  return sites_before_other_images + sites_before_within_image;
}

static int SiteIDNumberFromGrid(int x, int y, int z, int xResolution, int yResolution, int zResolution) {
  int sites_per_image = (int) (ceil(GSV_WIDTH / (double) scale_factor) * ceil(GSV_HEIGHT / (double) scale_factor));
  int sites_from_images = sites_per_image * number_images;
  int sites_before_within_grid = x + (xResolution) * y + (xResolution + yResolution) * z;
  return sites_from_images + sites_before_within_grid;
}

// generate grids for object files
static int GenerateGraph(void) {
  
  int converted_run = ConvertIndividualRun(scene, run_number);

  // find bounding box
  R3Box boundingBox(1, 1, 1, 0, 0, 0);
  boundingBox.Union(scene->Run(converted_run)->BBox());
  double grid_resolution = 0.67;

  // create grid
  int xLength = (int) ceil(ceil(boundingBox.XLength()) * grid_resolution);
  int yLength = (int) ceil(ceil(boundingBox.YLength()) * grid_resolution);
  int zLength = (int) ceil(ceil(boundingBox.ZLength()) * grid_resolution);

  R3Grid *grid = new R3Grid(xLength, yLength, zLength, boundingBox);
  
  // find the number of sites
  number_images = scene->Run(converted_run)->NImages();
  int sites_from_images = (int) (number_images * ceil(GSV_WIDTH / (double) scale_factor) * ceil(GSV_HEIGHT / (double) scale_factor));
  int sites_from_grids = grid->NEntries();
  int num_sites = sites_from_images + sites_from_grids;

  // create GCoptimization
  GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_sites, 2);

  
  for (int x = 0; x < xLength; x++) {
    for (int y = 0; y < yLength; y++) {
      for (int z = 0; z < zLength; z++) {
	int site_id = SiteIDNumberFromGrid(x, y, z, xLength, yLength, zLength);
	gc->setDataCost(site_id, 0, 0);
	gc->setDataCost(site_id, 1, 0);
      }
    }
  }
    
  int number_of_images_seen = 0;

  // read all of the files in the directory
  GSVRun *run = scene->Run(converted_run);
  for (int is = 0; is < run->NSegments(); is++) {
    GSVSegment *segment = run->Segment(is);
    for (int ip = 0; ip < segment->NPanoramas(); ip++) {
      GSVPanorama *panorama = segment->Panorama(ip);
      for (int ii = 0; ii < panorama->NImages(); ii++) {
	GSVImage *image = panorama->Image(ii);
	
	// find the file
	char vocFilename[128];
	snprintf(vocFilename, 127, "%s\%02d_%06d_%02d_UndistortedImage.txt", voc_directory, is, ip, ii);
	
	ifstream voc_probabilities;
	voc_probabilities.open(vocFilename);;
	if (voc_probabilities.is_open()) {
	  string line;
	  getline(voc_probabilities, line);
	  
	  // find scale factor and create array for all probabilities
	  int found_scale = atoi(line.c_str());
	  if (scale_factor != found_scale) {
	    fprintf(stderr, "Global scale factor (%d) does not match found scale (%d)", scale_factor, found_scale);
	    return 0;
	  }
	  
	  
	  double probability_array[(int) ceil(GSV_WIDTH / (double) scale_factor)][(int) ceil(GSV_HEIGHT / (double) scale_factor)];
	  
	  // for debugging purposes
	  int line_number = 1;
	  // read in all probabilities
	  while (getline(voc_probabilities, line)) {
	    vector<string> probabilities = splitString(line, ',');
	    
	    // make sure each line is formatted correctly
	    if (probabilities.size() != (unsigned int) ceil((unsigned int) GSV_HEIGHT / (double) scale_factor)) {
	      fprintf(stderr, "Incorrectly formatted line %d\n", line_number);
	      fprintf(stderr, "Found %u probabilities, expecting %u\n", probabilities.size(), (unsigned int) GSV_HEIGHT / scale_factor);
	      return 0;
	    }

	    // go through each probability
	    for (unsigned int i = 0; i < probabilities.size(); i++) {
	      probability_array[line_number - 1][i] = atof(probabilities[i].c_str());
	    }
	    line_number++;
	  }
	  
	  printf("Read %s\n", vocFilename);
	
	  // iterate through probability array
	  for (int y = 0; y < (int) (ceil(GSV_HEIGHT / (double) scale_factor)); y++) {
	    for (int x = 0; x < (int) (ceil(GSV_WIDTH / (double) scale_factor)); x++) {
	      int site_id = SiteIDNumberFromImage(number_of_images_seen, x, y);
	      
	      // set the data cost for this pixel
	      double weight = 50.0;
	      gc->setDataCost(site_id, 0, (weight * probability_array[x][y]));
	      gc->setDataCost(site_id, 1, (weight - weight * probability_array[x][y]));
	      
	      if (probability_array[x][y] > 0.5) {
		gc->setLabel(site_id, 1);
	      }

	      // create neighbors between pixels
	      if (x != 0) {
		int x_low = x - 1;
		int neighbor_id = SiteIDNumberFromImage(number_of_images_seen, x_low, y);
		gc->setNeighbors(site_id, neighbor_id);
	      }
	      if (y != 0) {
		int y_low = y - 1;
		int neighbor_id = SiteIDNumberFromImage(number_of_images_seen, x, y_low);
		gc->setNeighbors(site_id, neighbor_id);
	      }
	      
	      // create neighbors with appropriate grid cells
	      double xpoint = x * scale_factor;
	      double ypoint = y * scale_factor;
	      
	      R3Ray ray = image->RayThroughUndistortedPosition(R2Point((double) xpoint, (double) ypoint));
	      R3Point start = ray.Point(0);
	      R3Point end;
	      R3Intersects(ray, boundingBox, &end);
	      
	      R3Point gridStart = grid->GridPosition(start);
	      R3Point gridEnd = grid->GridPosition(end);
	      
	      double *pre1 = (double *) gridStart.Coords();
	      double *pre2 = (double *) gridEnd.Coords();
	      int p1[3] = { (int) (pre1[0] + 0.5), (int) (pre1[1] + 0.5), (int) (pre1[2] + 0.5) };
	      int p2[3] = { (int) (pre2[0] + 0.5), (int) (pre2[1] + 0.5), (int) (pre2[2] + 0.5) };
	      
	      // Get some convenient variables
	      int d[3],p[3],dd[3],s[3];
	      for (int i = 0; i < 3; i++) {
		d[i]= p2[i] - p1[i];
		if(d[i]<0){
		  dd[i] = -d[i];
		  s[i] = -1;
		}
		else{
		  dd[i] = d[i];
		  s[i] = 1;
		}
		p[i] = p1[i];
	      }
	      
	      // Choose dimensions
	      int i1=0;
	      if(dd[1]>dd[i1]){i1=1;}
	      if(dd[2]>dd[i1]){i1=2;}
	      int i2=(i1+1)%3;
	      int i3=(i1+2)%3;
	      
	      // Check span extent
	      if(dd[i1]==0){
		// Span is a point - rasterize it
		if (((p[0] >= 0) && (p[0] < grid->XResolution())) &&
		    ((p[1] >= 0) && (p[1] < grid->YResolution())) &&
		    ((p[2] >= 0) && (p[2] < grid->ZResolution()))) {
		  
		  int grid_site = SiteIDNumberFromGrid(p[0], p[1], p[2], grid->XResolution(), grid->YResolution(), grid->ZResolution());
		  gc->setNeighbors(site_id, grid_site);
		}
	      }
	      else {
		// Step along span
		int off[3] = { 0, 0, 0 };
		for (int i = 0; i <= dd[i1]; i++) {
		  if (((p[0] >= 0) && (p[0] < grid->XResolution())) &&
		      ((p[1] >= 0) && (p[1] < grid->YResolution())) &&
		      ((p[2] >= 0) && (p[2] < grid->ZResolution()))) {
		    int grid_site = SiteIDNumberFromGrid(p[0], p[1], p[2], grid->XResolution(), grid->YResolution(), grid->ZResolution());
		    gc->setNeighbors(site_id, grid_site);
		  }
		  off[i2]+=dd[i2];
		  off[i3]+=dd[i3];
		  p[i1]+=s[i1];
		  p[i2]+=s[i2]*off[i2]/dd[i1];
		  p[i3]+=s[i3]*off[i3]/dd[i1];
		  off[i2]%=dd[i1];
		  off[i3]%=dd[i1];
		}
	      }
	      
	      
	    }
	  }
	  
	} else {
	  fprintf(stderr, "Failed to read %s...Aborting...\n", vocFilename);
	  return 0;
	}
	

	// close the file and increment the number of images seen
	voc_probabilities.close();
	number_of_images_seen++;
      }
    }
  }
  
  gc->setSmoothCost(0, 1, 10);
  printf("%lld\n", gc->compute_energy());
  printf("There\n");
  gc->expansion(100000);
  printf("%lld\n", gc->compute_energy());
  printf("There\n");


  // go through each grid cell and create label grid
  for (int x = 0; x < xLength; x++) {
    for (int y = 0; y < yLength; y++) {
      for (int z = 0; z < zLength; z++) {
	int site_id = SiteIDNumberFromGrid(x, y, z, xLength, yLength, zLength);
	grid->SetGridValue(x, y, z, gc->whatLabel(site_id));
      }
    }
  }

  WriteGrid("file.grd", grid);
  

  // delete the grid
  delete grid;



  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static void Usage(void) {
  fprintf(stderr, "[-v] <gsv_scene> <input_grid> <voc_directory> <output_filename> -scale_factor <scale> -run_number <run>\n");
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
  // parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene = *argv; }
      else if (!strcmp(*argv, "-voc_directory")) { argc--; argv++; voc_directory = *argv; }
      else if (!strcmp(*argv, "-input_grid")) { argc--; argv++; input_grid = *argv; }
      else if (!strcmp(*argv, "-output_filename")) { argc--; argv++; output_filename = *argv; }
      else if (!strcmp(*argv, "-run_number")) { argc--; argv++; run_number = atoi(*argv); }
      else if (!strcmp(*argv, "-usage")) { return 0; }
      else if (!strcmp(*argv, "-scale_factor")) { argc--; argv++; scale_factor = atoi(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!input_scene) { input_scene = *argv; }
      else if (!input_grid) { input_grid = *argv; }
      else if (!voc_directory) { voc_directory = *argv; }
      else if (!output_filename) { output_filename = *argv; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    argc--; argv++;
  }

  if (!input_scene) { fprintf(stderr, "Must specify input scene\n"); return 0; }
  if (!input_grid) { fprintf(stderr, "Must specify probability grid\n"); return 0; }
  if (!voc_directory) { fprintf(stderr, "Must specify voc probability directory\n"); return 0; }
  if (!output_filename) { fprintf(stderr, "Must specify an output filename\n"); return 0; }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0;}

  // read in the gsv scene
  if (!ReadScene()) { return 0; }

  // make sure the run number is valid
  if (run_number < 0 || run_number > scene->NRuns()) { 
    fprintf(stderr, "Run %d does not exist in this scene\n", run_number); 
    return 0; 
  }

  // Read in the grid to memory
  //if (!ReadGrid()) { return 0; }
  
  // generate the graph
  if (!GenerateGraph()) { return 0; }

  // return SUCCESS
  return 1;
}
