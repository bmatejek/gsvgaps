////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"
#include "RNBasics/RNBasics.h"
extern "C" {
#include <CSparse/cs.h>
};

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static int run_number = -1;
static int scale_factor = 0;
static int gsv_notated_run = -1;
static char *object = NULL;
static char *scene_name = NULL;
static char *output_grid_name = NULL;
static GSVScene *scene = NULL;
static R3Grid *optimized = NULL;

enum {
  GSV_WIDTH = 1936,
  GSV_HEIGHT = 2592,
};

// global variables
static int x_per_image = 0;
static int y_per_image = 0;
static int total_images = 0;
static int x_grid_cells = 0;
static int y_grid_cells = 0;
static int z_grid_cells = 0;

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////
// Create the optimal output grid
////////////////////////////////////////////////////////////////////////

static int VariableIDFromImage(int image_number, int x, int y) {
  int variables_from_previous_images = image_number * (x_per_image * y_per_image);
  int variables_from_this_image = (x_per_image) * y + x;
  return variables_from_previous_images + variables_from_this_image;
}

static int VariableIDFromGrid(int x, int y, int z) {
  int variables_from_images = total_images * (x_per_image * y_per_image);
  int variables_from_this_grid = (x_grid_cells * y_grid_cells) * z + (x_grid_cells) * y + x;
  return variables_from_images + variables_from_this_grid;
}

static int OptimizeGrid(void) {
  // keep track of the number of images
  int image_number = 0;

  GSVRun *run = scene->Run(gsv_notated_run);
  
  // find the total number of images
  total_images = run->NImages();

  // create grid from bounding box
  R3Box boundingBox(1, 1, 1, 0, 0, 0);
  boundingBox.Union(run->BBox());
  double grid_resolution = 0.25;
  
  x_grid_cells = (int) ceil(ceil(boundingBox.XLength()) * grid_resolution);
  y_grid_cells = (int) ceil(ceil(boundingBox.YLength()) * grid_resolution);
  z_grid_cells = (int) ceil(ceil(boundingBox.ZLength()) * grid_resolution);

  optimized = new R3Grid(x_grid_cells, y_grid_cells, z_grid_cells, boundingBox);

  // m is the number of equations
  int equations_between_pixels_and_voxels = (x_per_image * y_per_image) * total_images;
  int equations_between_pixels = (2 * x_per_image * y_per_image - y_per_image - x_per_image) * total_images;
  int equations_between_voxels = 3 * x_grid_cells * y_grid_cells * z_grid_cells - x_grid_cells * y_grid_cells - y_grid_cells * z_grid_cells - x_grid_cells * z_grid_cells;
  int equations_for_pixels = (x_per_image * y_per_image) * total_images;
  int m = equations_between_pixels_and_voxels + equations_between_pixels + equations_between_voxels + equations_for_pixels;
  
  // n is the number of variables
  int n = (x_per_image * y_per_image) * total_images + x_grid_cells * y_grid_cells * z_grid_cells;

  // maxnz is the max number of nonzero elements
  int max_diagonal_grid_path = (int) ceil(sqrt(x_grid_cells * x_grid_cells + y_grid_cells * y_grid_cells + z_grid_cells * z_grid_cells));
  int nonzero_from_between_pixels_and_voxels = equations_between_pixels_and_voxels * max_diagonal_grid_path;
  int nonzero_from_between_pixels = 2 * equations_between_pixels;
  int nonzero_from_between_voxels = 2 * equations_between_voxels;
  int nonzero_for_pixels = equations_for_pixels;
  int maxnz = nonzero_from_between_pixels_and_voxels + nonzero_from_between_pixels + nonzero_from_between_voxels + nonzero_for_pixels;

  // create both needed matricies
  cs *M = cs_spalloc(m, n, maxnz, 1, 1);
  if (!M) {
    fprintf(stderr, "Failed to allocate matrix\n");
    return 0;
  }
  double *b = new double[m];
  for (int i = 0; i < m; i++) {
    b[i] = 0.0;
  }


  // keep track of how many equations have been input
  int number_equations_input = 0;

  // add in grid smoothing terms
  for (int x = 0; x < x_grid_cells; x++) {
    for (int y = 0; y < y_grid_cells; y++) {
      for (int z = 0; z < z_grid_cells; z++) {
	
	int this_variable_id = VariableIDFromGrid(x, y, z);
	if (this_variable_id > n) {
	  fprintf(stderr, "14: Error\n");
	  return 0;
	}
	
	// add all 3 smoothing terms
	if (x != x_grid_cells - 1) {
	  int neighbor_variable_id = VariableIDFromGrid(x + 1, y, z);
	  
	  // ADD IN SMOOTHING TERM
	  cs_entry(M, number_equations_input, this_variable_id, 1);
	  if (this_variable_id > n) {
	    fprintf(stderr, "13: Error\n");
	    return 0;
	  }
	  cs_entry(M, number_equations_input, neighbor_variable_id, -1);
	  if (neighbor_variable_id > n) {
	    fprintf(stderr, "12: Error\n");
	    return 0;
	  }
	  number_equations_input++;
	}
	
	if (y != y_grid_cells - 1) {
	  int neighbor_variable_id = VariableIDFromGrid(x, y + 1, z);

	  // ADD IN SMOOTHING TERM
	  cs_entry(M, number_equations_input, this_variable_id, 1);
	  if (this_variable_id > n) {
	    fprintf(stderr, "11: Error\n");
	    return 0;
	  }
	  cs_entry(M, number_equations_input, neighbor_variable_id, -1);
	  if (neighbor_variable_id > n) {
	    fprintf(stderr, "10: Error\n");
	    return 0;
	  }
	  number_equations_input++;
	}

	if (z != z_grid_cells - 1) {
	  int neighbor_variable_id = VariableIDFromGrid(x, y, z + 1);

	  // ADD IN SMOOTHING TERM
	  cs_entry(M, number_equations_input, this_variable_id, 1);
	  if (this_variable_id > n) {
	    fprintf(stderr, "9: Error\n");
	    return 0;
	  }
	  cs_entry(M, number_equations_input, neighbor_variable_id, -1);
	  if (neighbor_variable_id > n) {
	    fprintf(stderr, "8: Error\n");
	    return 0;
	  }
	  number_equations_input++;
	}
      }
    }
  }

  // iterate through all images
  for (int is = 0; is < run->NSegments(); is++) {
    GSVSegment *segment = run->Segment(is);
    for (int ip = 0; ip < segment->NPanoramas(); ip++) {
      GSVPanorama *panorama = segment->Panorama(ip);
      for (int ii = 0; ii < panorama->NImages(); ii++) {
	GSVImage *image = panorama->Image(ii);

	// add all smoothing terms for the images
	for (int x = 0; x < x_per_image; x++) {
	  for (int y = 0; y < y_per_image; y++) {
	    int this_variable_id = VariableIDFromImage(image_number, x, y);
	    
	    if (x != x_per_image - 1) {
	      int neighbor_variable_id = VariableIDFromImage(image_number, x + 1, y);

	      // ADD IN SMOOTHING TERM
	      cs_entry(M, number_equations_input, this_variable_id, 1);
	      if (this_variable_id > n) {
		fprintf(stderr, "7: Error\n");
		return 0;
	      }
	      cs_entry(M, number_equations_input, neighbor_variable_id, -1);
	      if (neighbor_variable_id > n) {
		fprintf(stderr, "6: Error\n");
		return 0;
	      }
	      number_equations_input++;
	    }

	    if (y != y_per_image - 1) {
	      int neighbor_variable_id = VariableIDFromImage(image_number, x, y + 1);

	      // ADD IN SMOOTHING TERM
	      cs_entry(M, number_equations_input, this_variable_id, 1);
	      if (this_variable_id > n) {
		fprintf(stderr, "5: Error\n");
		return 0;
	      }
	      cs_entry(M, number_equations_input, neighbor_variable_id, -1);
	      if (neighbor_variable_id > n) {
		fprintf(stderr, "4: Error\n");
		return 0;
	      }
	      number_equations_input++;
	    }
	  }
	}

	char probFilename[128];
	snprintf(probFilename, 127, "voc_probabilities/%s/%02d/%02d/%02d_%06d_%02d_UndistortedImage.txt", object, run_number, scale_factor, is, ip, ii);

	ifstream probabilities(probFilename);
	if (probabilities.is_open()) {
	  string line;
	  getline(probabilities, line);

	  // check scale factor
	  int scale = atoi(line.c_str());
	  if (scale != scale_factor) {
	    fprintf(stderr, "For %s, scale does not match indicated scale factor\n", probFilename);
	  }

	  // each line represents an x value
	  int x = 0;
	  while (getline(probabilities, line)) {
	    vector<string> scores = splitString(line, ',');
	    // each score corresponds to a different 
	    for (unsigned int y = 0; y < scores.size(); y++) {
	      int this_variable_id = VariableIDFromImage(image_number, x, y);

	      double probability = atof(scores[y].c_str());

	      // ADD IN DATA TERM FOR PIXEL
	      cs_entry(M, number_equations_input, this_variable_id, 1);
	      if (this_variable_id > n) {
		fprintf(stderr, "3: Error\n");
		return 0;
	      }
	      b[number_equations_input] = probability;
	      number_equations_input++;
	      
	      double xpoint = x * scale_factor;
	      double ypoint = y * scale_factor;

	      // find world points
	      R3Ray ray = image->RayThroughUndistortedPosition(R2Point(xpoint, ypoint));
	      R3Point start = ray.Point(0);
	      R3Point end;
	      R3Intersects(ray, boundingBox, &end);
	      
	      // convert world points to grid points
	      R3Point gridStart = optimized->GridPosition(start);
	      R3Point gridEnd = optimized->GridPosition(end);
	      
	      // add all data terms by iterating through points


	      // CODE MODIFIED FROM R3GRID
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
		if (((p[0] >= 0) && (p[0] < x_grid_cells)) &&
		    ((p[1] >= 0) && (p[1] < y_grid_cells)) &&
		    ((p[2] >= 0) && (p[2] < z_grid_cells))) {
		  
		  int variable_grid_id = VariableIDFromGrid(p[0], p[1], p[2]);
		  
		  // ADD IN DATA TERM
		  cs_entry(M, number_equations_input, variable_grid_id, 1);
		  if (variable_grid_id > n) {
		    fprintf(stderr, "2: Error\n");
		    return 0;
		  }
		}
	      }
	      else {
		// Step along span
		int off[3] = { 0, 0, 0 };
		for (int i = 0; i <= dd[i1]; i++) {
		  if (((p[0] >= 0) && (p[0] < x_grid_cells)) &&
		      ((p[1] >= 0) && (p[1] < y_grid_cells)) &&
		      ((p[2] >= 0) && (p[2] < z_grid_cells))) {
		    int variable_grid_id = VariableIDFromGrid(p[0], p[1], p[2]);

		    // ADD IN DATA TERM
		    cs_entry(M, number_equations_input, variable_grid_id, 1);
		    if (variable_grid_id > n) {
		      fprintf(stderr, "1: Error\n");
		      return 0;
		    }
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
	      b[number_equations_input] = probability;
	      number_equations_input++;
	    }
	    x++;
	  }
	  
	  printf("Read %s\n", probFilename);

	}
	else {
	  fprintf(stderr, "Unable to open %s\n", probFilename);
	  return 0;
	}
      }
    }
  }

  printf("%d %d\n", m, number_equations_input);

  // SOLVE LINEAR PROGRAM
  double *x = new double[n];
  if (!x) {
    fprintf(stderr, "Failed to allocate memory for results vector\n");
    return 0;
  }
  
  // set up aT * a * x = aT * b
  cs *A = cs_compress(M);
  cs_spfree(M);
  
  cs *AT = cs_transpose(A, 1);
  cs *ATA = cs_multiply(AT, A);
    
  cs_gaxpy(AT, b, x);
  
  if (!cs_cholsol(1, ATA, x)) {
    fprintf(stderr, "Error in cs_cholsol\n");
    return 0;
  }
  printf("x\n");
  for (int i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }

  // CREATE GRID BASED ON RESULTS

  // free all data
  cs_spfree(A);
  cs_spfree(AT);
  cs_spfree(ATA);
  delete [] b;
  delete [] x;

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static void Usage(void) {
  fprintf(stderr, "<gsv_scene> <output_grid> [-v] -object <object> -scale_factor <scale> -run_number <run_number>\n");
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
  // parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-object")) { argc--; argv++; object = *argv; }
      else if (!strcmp(*argv, "-scale_factor")) { argc--; argv++; scale_factor = atoi(*argv); }
      else if (!strcmp(*argv, "-run_number")) { argc--; argv++; run_number = atoi(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!scene_name) { scene_name = *argv; }
      else if (!output_grid_name) { output_grid_name = *argv; }
    }
    argc--; argv++;
  }
  
  // make sure all requirements are there 
  if (!scene_name) { fprintf(stderr, "Need to provide a gsv scene\n"); return 0; }
  if (!output_grid_name) { fprintf(stderr, "Need to provide an output grid name\n"); return 0; }
  if (!object) { fprintf(stderr, "Need to provide which object to analyze\n"); return 0; }
  if (!scale_factor) { fprintf(stderr, "Need to provide a scale factor\n"); return 0; }
  if (run_number == -1) { fprintf(stderr, "Need to provide a run number\n"); return 0; }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0;}

  x_per_image = (int) ceil(GSV_WIDTH / (double) scale_factor);
  y_per_image = (int) ceil(GSV_HEIGHT / (double) scale_factor);

  // read in the input scene
  scene = ReadScene(scene_name);
  if (!scene) { fprintf(stderr, "Failed to read input scene: %s\n", scene_name); return 0; }

  gsv_notated_run = ConvertIndividualRun(scene, run_number);
  
  if (!OptimizeGrid()) { fprintf(stderr, "Failed to create an optimized grid\n"); return 0; }

  //WriteGrid(output_grid_name, optimized);
  
  // return SUCCESS
  return 1;
}
