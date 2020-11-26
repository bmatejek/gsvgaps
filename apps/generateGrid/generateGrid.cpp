////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "generateGrid.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static char *input_scene_name = NULL;
static char *objects_of_interest_filename = NULL;
static char *output_grid_directory = NULL;
vector<int> run_input = vector<int>();

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

// write grid
static void WriteGrid(const char *filename, R3Grid *grid) {
  RNTime start_time;
  start_time.Read();
  
  if (!grid->WriteGridFile(filename)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return;
  }

  if (print_verbose) {
    printf("Write grid to %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
}

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

////////////////////////////////////////////////////////////////////////
// Helper Function
////////////////////////////////////////////////////////////////////////

// convert the point from the grid->BBox() into absolute grid coordinates
static R3Point AdjustToGrid(R3Point point, R3Grid *grid) {
  R3Box bbox = grid->WorldBox();
  double xDist = point.X() - bbox.XMin();
  double yDist = point.Y() - bbox.YMin();
  double zDist = point.Z() - bbox.ZMin();
  
  double xProportion = xDist / bbox.XLength();
  double yProportion = yDist / bbox.YLength();
  double zProportion = zDist / bbox.ZLength();
  
  double adjustedX = xProportion * grid->XResolution();
  double adjustedY = yProportion * grid->YResolution();
  double adjustedZ = zProportion * grid->ZResolution();

  return R3Point(adjustedX, adjustedY, adjustedZ);
}

////////////////////////////////////////////////////////////////////////
// ObjectOfInterest class functions
////////////////////////////////////////////////////////////////////////

// constructors
ObjectOfInterest ::
ObjectOfInterest(char *object_name, int index) :
  object_name(object_name),
  height(0.0),
  index(index)
{
}

// accessors/manipulators
void ObjectOfInterest ::
SetHeight(double height) {
  this->height = height;
}

const char *ObjectOfInterest ::
ObjectName(void) const {
  return object_name;
}

const double ObjectOfInterest ::
Height(void) const {
  return height;
}

const int ObjectOfInterest ::
Index(void) const {
  return index;
}


////////////////////////////////////////////////////////////////////////
// GenerateGrid class functions
////////////////////////////////////////////////////////////////////////

// constructors
GenerateGrid ::
GenerateGrid(GSVScene *scene, GSVPoseOptimization *optimization, bool *show_runs, vector<ObjectOfInterest *> objects) :
scene(scene),
grid(vector<R3Grid *>()),
optimization(optimization),
show_runs(show_runs),
objects_of_interest(objects)
{

  // create bounding box over visible runs
  R3Box bbox = R3Box(1, 1, 1, 0, 0, 0);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_runs[ir]) {
      bbox.Union(scene->Run(ir)->BBox());
    }
  }

  double xLength = bbox.XLength();
  double yLength = bbox.YLength();
  double zLength = bbox.ZLength();

  int gridXLength, gridYLength, gridZLength;
  // create grid from scene bounding boxes
  /*  if (xLength > 1000 && xLength < yLength) {
    gridXLength = 1000;
    gridYLength = (int) ceil(yLength * 1000 / xLength);
    gridZLength = (int) ceil(zLength * 1000 / xLength);
  }
  else if (yLength > 1000) {
    gridXLength = (int) ceil(xLength * 1000 / yLength);
    gridYLength = 1000;
    gridZLength = (int) ceil(zLength * 1000 / yLength);
  }
  else {*/
  gridXLength = (int) ceil(xLength);
  gridYLength = (int) ceil(yLength);
  gridZLength = (int) ceil(zLength);
    //}
  
  R3Grid *imageGrid = new R3Grid(gridXLength, gridYLength, gridZLength, bbox);
  if (!imageGrid) {
    fprintf(stderr, "Failed to allocate memory for image grid\n");
    return;
  }

  // write grid file for image locations
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_runs[ir]) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
	GSVSegment *segment = run->Segment(is);
	for (int ip = 0; ip < segment->NPanoramas(); ip++) {
	  GSVPanorama *panorama = segment->Panorama(ip);
	  for (int ii = 0; ii < panorama->NImages(); ii++) {
	    GSVImage *image = panorama->Image(ii);
	    GSVPose pose = image->Pose();
	    R3Point viewpoint = pose.Viewpoint();
	    // find adjusted x, y, z position in world space
	    R3Point adjusted = AdjustToGrid(viewpoint, imageGrid);
	    
	    imageGrid->SetGridValue((int) floor(adjusted.X()), (int) floor(adjusted.Y()), (int) floor(adjusted.Z()), 1.0);
	  }
	}
      }
    }
  }

  // add grid to grid vector
  grid.push_back(imageGrid);

  // write the image grid
  string filename = string(output_grid_directory) + string("imageGrid.grd");
  WriteGrid(filename.c_str(), imageGrid);

  // write grid files for all objects
  for (unsigned int io = 0; io < objects.size(); io++) {
    ObjectOfInterest *object = objects[io];

    // create truth grid if found
    R3Grid *truthGrid = new R3Grid(gridXLength, gridYLength, gridZLength, bbox);
    if (!truthGrid) {
      fprintf(stderr, "Failed to allocate memory for truth grid: %s\n", object->ObjectName());
      return;
    }
    
    // get all truth points for this object
    ifstream truth ((string(object->ObjectName()) + string("/truth_data.txt")).c_str());
    if (truth.is_open()) {
      string line;
      while (getline(truth, line)) {
	vector<string> coordinates = splitString(line, ',');
	
	// create new points
	double x = atof(coordinates[0].c_str());
	double y = atof(coordinates[1].c_str());
	double z = atof(coordinates[2].c_str());

	R3Point objectPoint = R3Point(x, y, z);
	R3Point adjusted = AdjustToGrid(objectPoint, truthGrid);

	truthGrid->SetGridValue((int) floor(adjusted.X()), (int) floor(adjusted.Y()), (int) floor(adjusted.Z()), 1.0);
      }
      truth.close();
    }
    else {
      fprintf(stderr, "Truth file not found for %s ... Attempting to continue program by ignoring truth points.\n", object->ObjectName());
    }
    
    string truthFilename = string(output_grid_directory) + string(object->ObjectName()) + string("_truth.grd");
    WriteGrid(truthFilename.c_str(), truthGrid);
    

    // create grid file for discoverd objects
    R3Grid *objectGrid = new R3Grid(gridXLength, gridYLength, gridZLength, bbox);
    if (!objectGrid) {
      fprintf(stderr, "Failed to allocate memory for object grid: %s\n", object->ObjectName());
      return;
    }


    for (int ir = 0; ir < scene->NRuns(); ir++) {
      if (show_runs[ir]) {
	GSVRun *run = scene->Run(ir);
	int run_number = IndividualRunToHumanReadable(scene, ir);
	for (int is = 0; is < run->NSegments(); is++) {
	  GSVSegment *segment = run->Segment(is);
	  for (int ip = 0; ip < segment->NPanoramas(); ip++) {
	    GSVPanorama *panorama = segment->Panorama(ip);
	    // subtract one and ignore the sky image
	    for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
	      GSVImage *image = panorama->Image(ii);

	      // open up corresponding text file to find stop signs
	      char file[128];
	      snprintf(file, 128, "%s/%02d/%02d_%06d_%02d_UndistortedImage.txt", object->ObjectName(), run_number, is, ip, ii);
	      ifstream boxes(file);
	      if (boxes.is_open()) {
		// get parameters of each box
		string line;
		while (getline(boxes, line)) {
		  vector<string> parameters = splitString(line, ',');
		  
		  // parse parameters
		  double xmin = atof(parameters[0].c_str());
		  double ymin = atof(parameters[1].c_str());
		  double xmax = atof(parameters[2].c_str());
		  double ymax = atof(parameters[3].c_str());
		  double score = atof(parameters[5].c_str());
		  
		  VOCBox box = VOCBox(xmin, ymin, xmax, ymax, score, image, object);
		  
		  R3Point objectPoint = box.GlobalPoint();
		  // find adjusted x, y, z position in world space
		  objectGrid->RasterizeWorldPoint(objectPoint, box.Score() + 2.0);
		}
		boxes.close();
	      }
	      else {
		fprintf(stderr, "Unable to open file: %s.  Attempting to continue program by ignoring image.\n", file);
	      }
	    }
	  }
	}
      }
    }
    // add grid to grid vector
    grid.push_back(objectGrid);

    // write this object's grid
    string filename = string(output_grid_directory) + string(object->ObjectName()) + string(".grd");
    WriteGrid(filename.c_str(), objectGrid);
  }
}

// accessors/manipulators
const GSVScene *GenerateGrid ::
Scene(void) const {
  return scene;
}

const vector<R3Grid *> GenerateGrid ::
Grid(void) const {
  return grid;
}

const GSVPoseOptimization *GenerateGrid ::
Optimization(void) const {
  return optimization;
}

const bool *GenerateGrid ::
ShowRuns(void) const {
  return show_runs;
}

const vector<ObjectOfInterest *> GenerateGrid ::
ObjectsOfInterest(void) const {
  return objects_of_interest;
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
      else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene_name = *argv; }
      else if (!strcmp(*argv, "-objects_of_interest")) { argc--; argv++; objects_of_interest_filename = *argv; }
      else if (!strcmp(*argv, "-output_grid_directory")) { argc--; argv++; output_grid_directory = *argv; }
      else if (!strcmp(*argv, "-show_runs")) {
	argc--; argv++;
	while (argc > 0 && (*argv)[0] != '-') {
	  run_input.push_back(atoi(*argv));
	  argc--; argv++;
	}
	// move back one argument since outer while loop will increment once more
	if (argc != 0)
	  argc++; argv--;
      }
      else if (!strcmp(*argv, "-usage")) { return 0; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!objects_of_interest_filename) { objects_of_interest_filename = *argv; }
      else if (!output_grid_directory) { output_grid_directory = *argv; }
      else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return 0; }
    }
    argc--; argv++;
  }
  
  // require input scene
  if (!input_scene_name) {
    fprintf(stderr, "Input Scene required.\n"); return 0;
  }

  // require input info file
  if (!objects_of_interest_filename) {
    fprintf(stderr, "Objects of Interest file required.\n"); return 0; 
  }

  // require output grid directory
  if (!output_grid_directory) {
    fprintf(stderr, "Output grid directory path required.\n"); return 0; 
  }

  // return SUCCESS
  return 1;
}

// parse the input detection data file
static vector<ObjectOfInterest *> ParseDetectionDataFile(void) {
  vector<ObjectOfInterest *> objects = vector<ObjectOfInterest *>();

  FILE *fp = fopen(objects_of_interest_filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed to read object of interest data file: %s\n", objects_of_interest_filename);
    return objects;
  }

  // read file
  char buff[256];

  // get all objects of interest
  while (fgets(buff, 256, fp)) {
    // ignore anything perceived as a comment
    if (buff[0] == '#') continue;

    // create objects of interest
    int bufferLength = strlen(buff);
    if (bufferLength <= 1) continue;
    char *object_name = (char *) malloc(sizeof(char) * bufferLength);
    strcpy(object_name, buff);
    object_name[bufferLength - 1] = '\0';
    ObjectOfInterest *object = new ObjectOfInterest(object_name, objects.size());

    // add to object of interst vector
    objects.push_back(object);
  }

  // close the file
  fclose(fp);

  // return vector of objects
  return objects;
}

// update the information of the object of interest (i.e. height, etc.)
static int UpdateObjectOfInterest(ObjectOfInterest *object) {
  // open up description file
  char fileLocation[128];
  
  strcpy(fileLocation, "./");
  strcat(fileLocation, object->ObjectName());
  strcat(fileLocation, "/description.txt");
  
  FILE *fp = fopen(fileLocation, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed to read object of interest %s description file\n", object->ObjectName());
    return 0;
  }

  // get all data regarding this particular object of interest 
  char buff[256];
  while (fscanf(fp, "%s", buff) == (unsigned int) 1) {
    if (!strcmp(buff, "HEIGHT:")) {
      double height;
      // make sure a valid height is given
      if (fscanf(fp, "%lf", &height) != (unsigned int) 1) {
	fprintf(stderr, "Invalid height in %s\n", fileLocation);
	return 0;
      }
      object->SetHeight(height);
    }
  }

  // close the file
  fclose(fp);

  // return SUCCESS
  return 1;
}


////////////////////////////////////////////////////////////////////////
// Local Helper Functions
////////////////////////////////////////////////////////////////////////

// convert runs from human ordered to gsv ordered
static void ConvertRuns(GSVScene *scene) {
  for (unsigned int ir = 0; ir < run_input.size(); ir++) {
    run_input[ir] = ConvertIndividualRun(scene, run_input[ir]);
  }
}


////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse function arguments
  if (!ParseArgs(argc, argv)) { Usage(); return 0;}
  
  // parse input data file to find all objects of interest - print error if no objects
  vector<ObjectOfInterest *> objects = ParseDetectionDataFile();
  if (!objects.size()) { fprintf(stderr, "At least one object of interest required!\n"); return 0; }

  // go through each object of interest and get required data
  for (unsigned int i = 0; i < objects.size(); i++) {
    if (!UpdateObjectOfInterest(objects[i])) return 0;
  }

  // read in the scene
  GSVScene *scene = ReadScene(input_scene_name);
  if (!scene) { fprintf(stderr, "Unable to create scene from %s\n", input_scene_name); return 0; }

  // convert runs from human ordered to gsv ordered
  ConvertRuns(scene);

  // create boolean array of run numbers to show
  bool *show_runs = (bool *) malloc(sizeof(bool) * scene->NRuns());
  if (!show_runs) { delete scene; return 0; }
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (run_input.size()) {
      show_runs[ir] = false;
    } 
    else {
      show_runs[ir] = true;
    }
  }

  for (unsigned int ui = 0; ui < run_input.size(); ui++) {
    // validate run number input
    if (run_input[ui] < 0 || run_input[ui] > scene->NRuns()) {
      fprintf(stderr, "Invalid run number: %d.  Possible run numbers range from %d to %d\n", run_input[ui], 0, scene->NRuns());
      delete show_runs; delete scene; return 0;
    }
    else {
      show_runs[run_input[ui]] = true;
    }
  }

  // create GSVPoseOptimization
  GSVPoseOptimization *optimization = new GSVPoseOptimization(scene);
  if (!optimization) { fprintf(stderr, "Failed to create GSVPoseOptimization\n"); delete show_runs; delete scene; return 0; }

  // generate a grid for the scene
  GenerateGrid *grid = new GenerateGrid(scene, optimization, show_runs, objects);
  if (!grid->Grid().size()) { fprintf(stderr, "Failed to create any grids\n"); delete optimization; delete show_runs; delete scene; return 0; }
  
  // return SUCCESS
  return 1;
}
