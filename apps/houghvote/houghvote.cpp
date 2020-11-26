////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static int guess_distance = 1;
static int mask_non_maxima = 1;
static int path_grid = 0;
static int truth_grid = 0;
static int grid_resolution = 2;
static int blur_grid = 1;
static int ignore_cutoff = 0;
static char *input_scene_name = NULL;
static char *output_grid_directory = NULL;
static char *object_file = NULL;
static vector<int> run_numbers = vector<int>();
static vector<string> objects = vector<string>();

// program variables
static GSVScene *scene = NULL;
static GSVPoseOptimization *optimization = NULL;

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

// read the object file
static int ReadObjectFile(void) {
  RNTime start_time;
  start_time.Read();

  // read object_file
  ifstream objectFile(object_file);
  if (objectFile.is_open()) {
    string line;
    while (getline(objectFile, line)) {
      // ignore things perceived as comments
      if (line.find("#") != 0)
	objects.push_back(line);
      
    }
  }
  // return FAILURE if object_file is no read correctly
  else {
    fprintf(stderr, "Failed to read object file: %s\n", object_file);
    return 0;
  }
  
  // close the objectFile
  objectFile.close();

  // print statistics
  if (print_verbose) {
    printf("Read object file from %s ...\n", object_file);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Generate the Grids
////////////////////////////////////////////////////////////////////////

// create a grid from a gsv file
static R3Grid *gsv2grd(int xLength, int yLength, int zLength, R3Box boundingBox) {
  // create the gsv grid
  R3Grid *gsvGrid  = new R3Grid(xLength, yLength, zLength, boundingBox);
  if (!gsvGrid) { fprintf(stderr, "Failed to allocate memory for gsv grid\n"); return NULL; }
  for (unsigned int ir = 0; ir < run_numbers.size(); ir++) {
    GSVRun *run = scene->Run(run_numbers[ir]);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
	GSVScan *scan = segment->Scan(ia);
	// read in all scan points into the gsv scene
	if (!scan->ReadPoints()) {
	  fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), scan->RunIndex(), scan->SegmentIndex());
	  return NULL;
	}
	// read in all of the points and then release them
	for (int il = 0; il < scan->NScanlines(); il++) {
	  GSVScanline *scanline = scan->Scanline(il);
	  for (int ip = 0; ip < scanline->NPoints(); ip++) {
	    // get the point location and make that grid value equal to 1
	    R3Point worldPoint = scanline->PointPosition(ip);
	    R3Point gridPoint = gsvGrid->GridPosition(worldPoint);
	    gsvGrid->SetGridValue(gridPoint.X(), gridPoint.Y(), gridPoint.Z(), 1.0);
	  }
	}

	// release the points to save memory
	scan->ReleasePoints();
      }
    }
  }
  
  // return gsv grid
  return gsvGrid;
}

// create a grid for the car path
static R3Grid *path2grd(int xLength, int yLength, int zLength, R3Box boundingBox) {
  // for each image location, create a corresponding grid entry
  R3Grid *pathGrid = new R3Grid(xLength, yLength, zLength, boundingBox);
  if (!pathGrid) { fprintf(stderr, "Failed to allocate memory for path grid\n"); return NULL; }
  for (unsigned int ir = 0; ir < run_numbers.size(); ir++) {
    GSVRun *run = scene->Run(run_numbers[ir]);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
	GSVPanorama *panorama = segment->Panorama(ip);
	for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
	  GSVImage *image = panorama->Image(ii);
	  // get the image camera location and make that value on grid equal to 1
	  R3Point worldPoint = image->Pose().Viewpoint();
	  R3Point gridPoint = pathGrid->GridPosition(worldPoint);
	  pathGrid->SetGridValue(gridPoint.X(), gridPoint.Y(), gridPoint.Z(), 1.0);
	}
      }
    }
  }

  // return grid
  return pathGrid;
}

// create a truth grid for this object of interest
static R3Grid *truth2grd(int xLength, int yLength, int zLength, R3Box boundingBox, string object) {
  printf("Here\n");
  
  // find all truth locations
  vector<R3Point> truthPoints = vector<R3Point>();

  // create truth grid
  R3Grid *truthGrid = new R3Grid(xLength, yLength, zLength, boundingBox);
  if (!truthGrid) { fprintf(stderr, "Failed to allocate memory for %s truth grid\n", object.c_str()); return NULL; }
  
  // parse through truth file
  //string truthFilename = object + string("/truth_data.txt");
  string truthFilename = string("voc_predictions/traffic_light/truth_data.txt");
  printf("Here\n");
  printf("%s\n", truthFilename.c_str());
  ifstream truth("voc_predictions/traffic_light/truth_data.txt");
  printf("Why?");
  if (truth.is_open()) {
    string line;
    // get all truth poitns
    while (getline(truth, line)) {
      printf("line: %s", line.c_str());
      vector<string> coordinates = splitString(line, ',');

      // get world point coordinates
      double x = atof(coordinates[0].c_str());
      double y = atof(coordinates[1].c_str());
      double z = atof(coordinates[2].c_str());

      // get world point, convert to grid point, and add to truth grid
      R3Point worldPoint = R3Point(x, y, z);
      R3Point gridPoint = truthGrid->GridPosition(worldPoint);
      truthGrid->SetGridValue(gridPoint.X(), gridPoint.Y(), gridPoint.Z(), 1.0);
    }
  }
  // unable to open truth file
  else {
    printf("Why?\n");
    //fprintf(stderr, "Failed to open %s ... Unable to create truth grid for object %s\n", truthFilename.c_str(), object.c_str());
    return NULL;
  }
  
  // return grid
  return truthGrid;
}

// generate grids for object files
static int GenerateGrids(void) {
  // get bounding box
  R3Box boundingBox(1, 1, 1, 0, 0, 0);
  for (unsigned int i = 0; i < run_numbers.size(); i++) {
    boundingBox.Union(scene->Run(run_numbers[i])->BBox());
  }

  // create grid
  int xLength = ceil(boundingBox.XLength()) * grid_resolution;
  int yLength = ceil(boundingBox.YLength()) * grid_resolution;
  int zLength = ceil(boundingBox.ZLength()) * grid_resolution;
  
  // get all gsv points and convert gsv into grid
  R3Grid *gsvGrid = gsv2grd(xLength, yLength, zLength, boundingBox);
  if (!gsvGrid) { return 0; }

  // write gsv grid to file
  string gsvFilename = string(output_grid_directory) + string("gsvGrid.grd");
  WriteGrid(gsvFilename.c_str(), gsvGrid);
  
  // get the grid for the vehicle path
  if (path_grid) {
    R3Grid *pathGrid = path2grd(xLength, yLength, zLength, boundingBox);
    if (!pathGrid) { delete gsvGrid; return 0; }

    // write path grid to file
    string pathFilename = string(output_grid_directory) + string("pathGrid.grd");
    WriteGrid(pathFilename.c_str(), pathGrid);
    delete pathGrid;
  }
  
  // go through all objects of interest
  for (unsigned int io = 0; io < 1/*objects.size()*/; io++) {
    // get the object name
    string objectName = string("traffic_light"); // objects[io];

    // get truth data
    if (truth_grid) {
      R3Grid *truthGrid = truth2grd(xLength, yLength, zLength, boundingBox, objectName); 
      if (truthGrid) {
	// write truth grid to file
	//string truthFilename = string(output_grid_directory) + string("truth_") + objectName + string(".grd");
	//WriteGrid(truthFilename.c_str(), truthGrid);
	WriteGrid("truth_traffic_light.grd", truthGrid);
	delete truthGrid;
      }
    }
    
    // get voc data
    R3Grid *objGrid = new R3Grid(xLength, yLength, zLength, boundingBox);
    if (!objGrid) { fprintf(stderr, "Failed to allocate memory for grid\n"); delete gsvGrid; return 0; }
    
    // get the object parameters
    double objectHeight = -1.0, heightVariance = -1.0, vocVariance = -1.0, score_cutoff = -3.0;
    
    // open up the object description file
    char descriptionFile[128];
    snprintf(descriptionFile, 128, "probability_distributions/traffic_light/description.txt");
    
    ifstream description(descriptionFile);
    if (description.is_open()) {
      string line;
      while (getline(description, line)) {
	char *attribute = (char *) calloc(32, sizeof(char));
	double attribute_value;
	sscanf(line.c_str(), "%s%lf", attribute, &attribute_value);
	
	string attr = string(attribute);
	// get height attribute
	if (!attr.find("HEIGHT:")) {
	  objectHeight = attribute_value;
	}
	// get actual variance attribute
	else if (!attr.find("HEIGHT_VARIANCE:")) {
	  heightVariance = attribute_value;
	}
	// get voc variance
	else if (!attr.find("VOC_VARIANCE:")) {
	  vocVariance = attribute_value;
	}
	// get score cutoff
	else if (!attr.find("SCORE_CUTOFF:")) {
	  score_cutoff = attribute_value;
	}
	free(attribute);
      }
    }
    // give warning but continue program
    else {
      fprintf(stderr, "Failed to find: %s.  Ignoring this object...\n", descriptionFile);
      delete objGrid;
      continue;
    }

    // close the description file
    description.close();

    // make sure all description values are present
    bool pass = false;
    if (objectHeight < 0.0) {
      fprintf(stderr, "%s missing HEIGHT attribute!\n", descriptionFile);
      pass = true;
    }
    if (heightVariance < 0.0) {
      fprintf(stderr, "%s missing HEIGHT_VARIANCE attribute!\n", descriptionFile);
      pass = true;
    }
    if (vocVariance < 0.0) {
      fprintf(stderr, "%s missing VOC_VARIANCE attribute!\n", descriptionFile);
      pass = true;
    }
    if (!ignore_cutoff) {
      if (score_cutoff < -2.0) {
	fprintf(stderr, "%s missing SCORE_CUTOFF attribute (value must be greater than -2.0)!\n", descriptionFile);
	pass = true;
      }
    }
    else {
      score_cutoff = -0.5;
    }

    // skip object if value is missing
    if (pass) {
      fprintf(stderr, "Failed to create grid for %s.  Ignoring this object...\n", objectName.c_str());
      delete objGrid;
      continue;
    }

    // find lower and upper bounds
    double lowerBound = (1.0 - heightVariance) / (1.0 + vocVariance);
    double upperBound = (1.0 + heightVariance) / (1.0 - vocVariance);

    RNTime start_time;
    start_time.Read();

    // go through each run as input by user
    for (unsigned int ir = 0; ir < run_numbers.size(); ir++) {
      int run_number = IndividualRunToHumanReadable(scene, run_numbers[ir]);
      GSVRun *run = scene->Run(run_numbers[ir]);
      // go through each segment in run
      for (int is = 0; is < run->NSegments(); is++) {
	GSVSegment *segment = run->Segment(is);
	// go through each panorama in segment
	for (int ip = 0; ip < segment->NPanoramas(); ip++) {
	  GSVPanorama *panorama = segment->Panorama(ip);
	  // go through each image in panorama
	  for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
	    GSVImage *image = panorama->Image(ii);
	    
	    R3Vector towards = image->Pose().Towards();
	    R3Point viewpoint = image->Pose().Viewpoint();
	    
	    // open up corresponding file
	    char vocFilename[128];
	    snprintf(vocFilename, 127, "probability_distributions/traffic_light/%02d/scores/%02d_%06d_%02d_UndistortedImage.txt", run_number, is, ip, ii);
	    
	    char heightFilename[128];
	    snprintf(heightFilename, 127, "probability_distributions/traffic_light/%02d/heights/%02d_%06d_%02d_UndistortedImage.txt", run_number, is, ip, ii);

	    ifstream heights(heightFilename);
	    if (!heights.is_open()) {
	      fprintf(stderr, "Failed to open %s\n", heightFilename);
	      return 0;
	    }
	    
	    // open up the individual voc filename
	    ifstream boxes(vocFilename);
	    if (boxes.is_open()) {
	      // get parameters for each box in the file
	      // go through each line
	      string line;
	      getline(boxes, line);

	      string lineHeight;
	      getline(heights, line);
	      
	      int scale_factor = atoi(line.c_str());
	      
	      // error for floating point numbers
	      double epsilon = 2; //10e-8;

	      int x = 0;
	      while (getline(boxes, line)) {
		getline(heights, lineHeight);
		vector<string> vocBox = splitString(line, ',');
		vector<string> heights = splitString(lineHeight, ',');

		
		int y = 0;
		for (unsigned int val = 0; val < vocBox.size(); val++) {
		  double score = atof(vocBox[val].c_str());
		  if (score > epsilon) {
		    // get height
		    double height = atof(heights[val].c_str());
		    
		    double image_plane_distance = EstimatedImagePlaneDistance(0, height, objectHeight, image);
		    double image_point_distance = EstimatedImagePointDistance(x, y, image_plane_distance, image);
		    
		    
		    R3Ray ray = image->RayThroughUndistortedPosition(R2Point((double) x, (double) y));
		    R3Point start = ray.Point(image_point_distance * lowerBound);
		    R3Point end = ray.Point(image_point_distance * upperBound);
		    //R3Point start = ray.Point(1.0);
		    //R3Point end;
		    R3Intersects(ray, boundingBox, &end);
		    objGrid->RasterizeWorldSpanHough(start, end, score - 2, viewpoint, towards);
		  }
		  y += scale_factor;
		}
		x += scale_factor;
	      }
	      
	      // add tmp to the objGrid
	      printf("Read %s\n", vocFilename);
	    }
	    // attempt to continue the program
	    else {
	      fprintf(stderr, "Failed to open file: %s....Attempting to recover....\n", vocFilename);
	    }
	  }
	}
      }
    }

    if (print_verbose) {
      printf("Read all files in %.2f seconds...\n", start_time.Elapsed());
    }
    
    // mask non maxima if desired
    if (mask_non_maxima) {
      objGrid->MaskNonMaxima();
    }
    
    // copy gsvGrid into mergeGrid to merge and delete
    R3Grid *mergeGrid = new R3Grid(*gsvGrid);
    *mergeGrid *= *objGrid;
    // different filename depending on whether or not the cutoff is ignored
    string mergeFilename;
    if (ignore_cutoff)
      mergeFilename = string(output_grid_directory) + string("no_cutoff_merged_") + string("traffic_light") string(".grd");
    else
      mergeFilename = string(output_grid_directory) + string("merged_") + string("traffic_light") + string(".grd");
    WriteGrid(mergeFilename.c_str(), mergeGrid);
    delete mergeGrid;

    // blur the grid if desired
    if (blur_grid)
      objGrid->Blur();  

    // write object grid
    string objectFilename = string(output_grid_directory) + string("traffic_light")  + string(".grd");
    WriteGrid(objectFilename.c_str(), objGrid);
    
    // delete the grid
    delete objGrid;
  }
  
  // delete gsvGrid
  delete gsvGrid;

  }
  // return SUCCESS
  return 1;
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
  // parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-dont_guess_distance")) { guess_distance = 0; }
      else if (!strcmp(*argv, "-dont_mask_non_maxima")) { mask_non_maxima = 0; }      
      else if (!strcmp(*argv, "-dont_blur")) { blur_grid = 0; }
      else if (!strcmp(*argv, "-draw_path")) { path_grid = 1; }
      else if (!strcmp(*argv, "-draw_truth")) { truth_grid = 1; }
      else if (!strcmp(*argv, "-ignore_cutoff")) { ignore_cutoff = 1; }
      else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene_name = *argv; }
      else if (!strcmp(*argv, "-output_grid_directory")) { argc--; argv++; output_grid_directory = *argv; }
      else if (!strcmp(*argv, "-grid_resolution")) { argc--; argv++; grid_resolution = atoi(*argv); }
      else if (!strcmp(*argv, "-usage")) { return 0; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_grid_directory) output_grid_directory = *argv;
      else if (!object_file) object_file = *argv;
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    argc--; argv++;
  }
  
  // require input scene
  if (!input_scene_name) { fprintf(stderr, "Input Scene required.\n"); return 0; }

  // have default output_grid_directory by the current working directory
  if (!output_grid_directory) { 
    fprintf(stderr, "Warning: output directory defaults to current working directory\n"); 
    output_grid_directory = (char *) "./";
  }

  // warn if no object file list
  if (!object_file) { 
    if (path_grid) 
      fprintf(stderr, "Warning: no object file specified, only creating %sgsvGrid.grd and %spathGrid.grd\n", output_grid_directory, output_grid_directory);
    else
      fprintf(stderr, "Warning: no object file specified, only creating %sgsvGrid.grd\n", output_grid_directory); 
  }

  // return SUCCESS
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0;}
  
  // read scene
  scene = ReadScene(input_scene_name);
  if (!scene) { fprintf(stderr, "Unable to read scene: %s\n", input_scene_name); return 0; }

  // create gsv optimization
  optimization = new GSVPoseOptimization(scene);
  if (!optimization) { fprintf(stderr, "Unable to create GSVPoseOptimization\n"); delete scene; return 0; }

  // parse through object file to find all objects of interest
  if (object_file && !ReadObjectFile()) { fprintf(stderr, "Failed to read object file %s\n", object_file); delete scene; delete optimization; return 0; }

  // convert all run numbers
  run_numbers.push_back(ConvertIndividualRun(scene, 10));

  // generate the desired grids
  if (!GenerateGrids()) { fprintf(stderr, "Failed to produce all desired grids\n"); delete scene; delete optimization; return 0; }

  // delete memory and exit
  delete optimization;
  delete scene;

  // return SUCCESS
  return 1;
}
