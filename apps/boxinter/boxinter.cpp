////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// program variables
static int print_verbose = 0;
static char *input_scene_name = NULL;
static char *objects_file_name = NULL;

// RCNN bounding box object
class RCNNBoundingBox {
	// constructor and access functions
public:
	RCNNBoundingBox(string object, R2Box bbox, double score, GSVImage *image);
	string Object(void);
	R2Box BBox(void);
	double Score(void);
	GSVImage *Image(void);

	// instance variables
private:
	string object;
	R2Box bbox;
	double score;
	GSVImage *image;
};

// constructor for RCNNBoundingBox
RCNNBoundingBox::RCNNBoundingBox(string object, R2Box bbox, double score, GSVImage *image) :
object(object),
bbox(bbox),
score(score),
image(image)
{
}

// return the object
string RCNNBoundingBox::Object(void) {
	return object;
}

// return the bounding box
R2Box RCNNBoundingBox::BBox(void) {
	return bbox;
}

// return the score of this RCNN bounding box
double RCNNBoundingBox::Score(void) {
	return score;
}

GSVImage *RCNNBoundingBox::Image(void) {
	return image;
}

// global variables
GSVScene *scene = NULL;								// GSVScene used through out the program
map<string, int> object_map;						// hash function to map object strings to ids
vector<RCNNBoundingBox *> ***obj_bboxes = NULL;		// bounding box array of all rcnn images
GSVImage **images = NULL;							// an array of all images

// useful constants
static int GSV_IMAGE_HEIGHT = 2592;


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *ReadScene(const char *filename) {
	// Start statistics
	RNTime start_time;
	start_time.Read();

	// Check if should read points (if .gsv, will be read as needed)
	RNBoolean read_points = (strstr(filename, ".gsv")) ? FALSE : TRUE;

	// Allocate google scene
	GSVScene *scene = new GSVScene();
	if (!scene) {
		fprintf(stderr, "Unable to allocate scene\n");
		return NULL;
	}

	// Read scene 
	if (!scene->ReadFile(filename, read_points)) {
		delete scene;
		return NULL;
	}

	// Print statistics
	if (print_verbose) {
		printf("Read scene from %s ...\n", filename);
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
		printf("  # Runs = %d\n", scene->NRuns());
		printf("  # Cameras = %d\n", scene->NCameras());
		printf("  # Lasers = %d\n", scene->NLasers());
		printf("  # Segments = %d\n", scene->NSegments());
		printf("  # Scans = %d\n", scene->NScans());
		printf("  # Tapestries = %d\n", scene->NTapestries());
		printf("  # Panoramas = %d\n", scene->NPanoramas());
		printf("  # Images = %d\n", scene->NImages());
		printf("  # Scanlines = %d\n", scene->NScanlines());
		fflush(stdout);
	}

	// Return scene
	return scene;
}

// split the line by delim and return the vector
static vector<string> splitString(string line, char delim = ',') {
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
// Find bounding box intersections
////////////////////////////////////////////////////////////////////////

static int ReadObjectFile(void) {
	// read the input file
	ifstream object_list;
	object_list.open(objects_file_name);

	// make sure the object is open
	if (!object_list.is_open()) { fprintf(stderr, "Unable to read object file at %s\n", objects_file_name); return 0; }

	// read in list of all objects 
	string line;
	while (getline(object_list, line)) {
		vector<string> object_data = splitString(line, ' ');

		// get id and object
		int id = atoi(object_data[0].c_str());
		string object = object_data[1];

		object_map.insert(pair<string, int>(object, id));
	}

	// return success
	return 1;
}

static int ReadAllImages(void) {
	images = new GSVImage *[scene->NImages()];
	if (!images) { fprintf(stderr, "Failed to allocate memory for image array\n"); return 0; }

	for (int ir = 0, in = 0; ir < scene->NRuns(); ++ir) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); ++is) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ++ii, ++in) {
					images[in] = panorama->Image(ii);
				}
			}
		}
	}

	// return success;
	return 1;
}

static int ReadBoundingBoxes(void) {
	// Start statistics
	RNTime start_time;
	start_time.Read();

	// create the array for all of the bounding boxes
	int nobjects = (int)object_map.size();
	int nimages = scene->NImages();

	// allocate memory for all vectors
	obj_bboxes = new vector<RCNNBoundingBox *> **[nimages];
	if (!obj_bboxes) { fprintf(stderr, "Failed to allocate memory for bounding boxes\n"); return 0; }
	for (int i = 0; i < nimages; i++) {
		obj_bboxes[i] = new vector<RCNNBoundingBox *> *[nobjects];
		if (!obj_bboxes[i]) { fprintf(stderr, "Failed to allocate memory for bounding boxes\n"); return 0; }
		for (int j = 0; j < nobjects; j++) {
			obj_bboxes[i][j] = new vector<RCNNBoundingBox *>();
		}
	}

	// go through all images
	int npredictions = 0;
	for (int ir = 0, in = 0; ir < scene->NRuns(); ++ir) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); ++is) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0, pn = 0; ip < segment->NPanoramas(); ++ip, ++pn) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ++ii, ++in) {
					GSVImage *image = panorama->Image(ii);

					// read the object file
					char object_file_name[4096];
					sprintf(object_file_name, "rcnn/%s/%02d_%06d_%02d_UndistortedImage.txt", run->Name(), is, pn, ii);
					ifstream object_file;
					object_file.open(object_file_name);

					if (!object_file.is_open()) { fprintf(stderr, "Unable to open %s\n", object_file_name); return 0; }

					// add in all of the bounding boxes for this file
					string line;
					while (getline(object_file, line)) {
						vector<string> object_info = splitString(line);

						// get all parameters for the RCNN prediction
						string object = object_info[0];
						double score = atof(object_info[1].c_str());
						double xmin = atof(object_info[2].c_str());
						double ymax = GSV_IMAGE_HEIGHT - atof(object_info[3].c_str());
						double xmax = atof(object_info[4].c_str());
						double ymin = GSV_IMAGE_HEIGHT - atof(object_info[5].c_str());

						RCNNBoundingBox *bbox = new RCNNBoundingBox(object, R2Box(xmin, ymin, xmax, ymax), score, image);

						// get the id from the object map
						int id = object_map[object];

						// add to the array of bounding boxes
						obj_bboxes[in][id]->push_back(bbox);
						npredictions++;
					}
				}
			}
		}
	}

	// Print statistics
	if (print_verbose) {
		printf("Read all RCNN object predictions ...\n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
		printf("  # Files read: %d\n", nimages);
		printf("  # Number of object classes: %d\n", nobjects);
		printf("  # Total number of RCNN predictions: %d\n", npredictions);
		fflush(stdout);
	}

	// return success
	return 1;
}

static void Intersection(RCNNBoundingBox *a, RCNNBoundingBox *b, ofstream *output) {
	// create four halfspaces for each bounding box
	GSVImage *ai = a->Image();
	GSVImage *bi = b->Image();
	R3Point aV = ai->Pose().Viewpoint();
	R3Point bV = bi->Pose().Viewpoint();

	double axmin = a->BBox().XMin(); double bxmin = b->BBox().XMin();
	double aymin = a->BBox().YMin(); double bymin = b->BBox().YMin();
	double axmax = a->BBox().XMax(); double bxmax = b->BBox().XMax();
	double aymax = a->BBox().YMax(); double bymax = b->BBox().YMax();

	R2Point a1 = R2Point(axmin, aymin); R2Point b1 = R2Point(bxmin, bymin);
	R2Point a2 = R2Point(axmin, aymax); R2Point b2 = R2Point(bxmin, bymax);
	R2Point a3 = R2Point(axmax, aymax); R2Point b3 = R2Point(bxmax, bymax);
	R2Point a4 = R2Point(axmax, aymin); R2Point b4 = R2Point(bxmax, bymin);

	R3Ray ra1 = ai->RayThroughUndistortedPosition(a1);
	R3Ray ra2 = ai->RayThroughUndistortedPosition(a2);
	R3Ray ra3 = ai->RayThroughUndistortedPosition(a3);
	R3Ray ra4 = ai->RayThroughUndistortedPosition(a4);

	R3Ray rb1 = bi->RayThroughUndistortedPosition(b1);
	R3Ray rb2 = bi->RayThroughUndistortedPosition(b2);
	R3Ray rb3 = bi->RayThroughUndistortedPosition(b3);
	R3Ray rb4 = bi->RayThroughUndistortedPosition(b4);

	R3Plane ppa1 = R3Plane(aV, ra1.Vector(), ra2.Vector());
	R3Plane ppa2 = R3Plane(aV, ra2.Vector(), ra3.Vector()); 
	R3Plane ppa3 = R3Plane(aV, ra3.Vector(), ra4.Vector());
	R3Plane ppa4 = R3Plane(aV, ra4.Vector(), ra1.Vector());

	R3Plane ppb1 = R3Plane(bV, rb1.Vector(), rb2.Vector());
	R3Plane ppb2 = R3Plane(bV, rb2.Vector(), rb3.Vector());
	R3Plane ppb3 = R3Plane(bV, rb3.Vector(), rb4.Vector());
	R3Plane ppb4 = R3Plane(bV, rb4.Vector(), rb1.Vector());

	*output << ppa1.A(); *output << ','; *output << ppa1.B(); *output << ','; *output << ppa1.C(); *output << ','; *output << ppa1.D();	*output << '\n';
	*output << ppa2.A(); *output << ','; *output << ppa2.B(); *output << ','; *output << ppa2.C(); *output << ','; *output << ppa2.D();	*output << '\n';
	*output << ppa3.A(); *output << ','; *output << ppa3.B(); *output << ','; *output << ppa3.C(); *output << ','; *output << ppa3.D();	*output << '\n';
	*output << ppa4.A(); *output << ','; *output << ppa4.B(); *output << ','; *output << ppa4.C(); *output << ','; *output << ppa4.D();	*output << '\n';
	
	*output << ppb1.A(); *output << ','; *output << ppb1.B(); *output << ','; *output << ppb1.C(); *output << ','; *output << ppb1.D();	*output << '\n';
	*output << ppb2.A(); *output << ','; *output << ppb2.B(); *output << ','; *output << ppb2.C(); *output << ','; *output << ppb2.D();	*output << '\n';
	*output << ppb3.A(); *output << ','; *output << ppb3.B(); *output << ','; *output << ppb3.C(); *output << ','; *output << ppb3.D();	*output << '\n';
	*output << ppb4.A(); *output << ','; *output << ppb4.B(); *output << ','; *output << ppb4.C(); *output << ','; *output << ppb4.D();	*output << '\n';
}

static int FindIntersections(void) {
	// Start statistics
	RNTime start_time;
	start_time.Read();

	ofstream *output = new ofstream();
	output->open("halfspaces.txt");

	if (!output->is_open()) { fprintf(stderr, "Failed to write to halfspaces.txt\n"); return 0; }

	// go through each object
	int nobjects = object_map.size();
	for (int io = 0; io < nobjects; io++) {
		// go through each unique pair of images
		for (int ii = 0; ii < scene->NImages(); ii++) {
			// get the bounding boxes for this image
			vector<RCNNBoundingBox *> *b = obj_bboxes[ii][io];
			for (int iip = 0; iip < ii; iip++) {
				// get the bounding boxes for this prime image
				vector<RCNNBoundingBox *> *bp = obj_bboxes[iip][io];

				// consider each pair of bounding boxes
				for (unsigned int ib = 0; ib < b->size(); ib++) {
					RCNNBoundingBox *bbox = (*b)[ib];
					if (bbox->Score() < 0.0) continue;
					for (unsigned int ibp = 0; ibp < bp->size(); ibp++) {
						RCNNBoundingBox *bboxp = (*bp)[ibp];
						if (bboxp->Score() < 0.0) continue;
						Intersection(bbox, bboxp, output);
					}
				}
			}
		}
	}

	// Print statistics
	if (print_verbose) {
		printf("Intersected all bounding boxes ...\n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
		fflush(stdout);
	}

	output->close();

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

// how should this program be used
static int Usage(void) {
	fprintf(stderr, "Usage: boxinter.cpp [-v] input_scene_name [-options ...]\n");
	return 0;
}

// argument parsing function
static int ParseArgs(int argc, char **argv) {
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-input_scene_name")) { input_scene_name = *argv; }
			else if (!strcmp(*argv, "-objects_file_name")) { objects_file_name = *argv; }
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else if (!objects_file_name) objects_file_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// make sure there is an input scene
	if (!input_scene_name) { fprintf(stderr, "Need input file\n"); return Usage(); }

	// set the default object file name
	if (!objects_file_name) { objects_file_name = "rcnn/objects.txt"; }

	// Return OK status 
	return 1;
}

// main method
int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// Read in scene
	scene = ReadScene(input_scene_name);
	if (!scene) { exit(-1); }

	// Read in all object file
	if (!ReadObjectFile()) { exit(-1); }

	// Read in all images
	if (!ReadAllImages()) { exit(-1); }

	// Read in all bounding boxes into an array
	if (!ReadBoundingBoxes()) { exit(-1); }

	// Aggregate the votes
	if (!FindIntersections()) { exit(-1); }

	// Return success 
	return 0;
}