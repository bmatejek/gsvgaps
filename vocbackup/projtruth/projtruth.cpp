// Source file for the google viewer program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include <vector>
#include <string>
#include <fstream>

using namespace std;



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments
static const char *input_scene_name = NULL;
static vector<char *> objects = vector<char *>();
static int print_verbose;

// global variables
static GSVScene *scene;

// constants
static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;

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
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage() {
	fprintf(stderr, "Usage: projtruth input_scene [objects]\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else { objects.push_back(*argv); }
			argv++; argc--;
		}
	}

	// Check all input
	if (!input_scene_name) return Usage();
	
	// make sure there are at least some objects
	if (objects.size() == 0) return Usage();

	// Return OK status 
	return 1;
}


////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// Read scene
	GSVScene *scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);

	// go through every object
	for (unsigned int io = 0; io < objects.size(); io++) {
		char *object = objects[io];

		vector<R3Point> true_locations = vector<R3Point>();

		// get every true bounding box
		for (int ir = 0; ir < scene->NRuns(); ir++) {
			GSVRun *run = scene->Run(ir);
			char truthFilename[4096];
			sprintf(truthFilename, "truth/%s/%s/truth.txt", object, run->Name());
		
			// open up truth file
			ifstream truth;
			truth.open(truthFilename);
			if (!truth.is_open()) {
				printf("Unable to open %s\n", truthFilename);
				continue;
			}

			// add in all the points
			string line;
			while (getline(truth, line)) {
				vector<string> coordinates = splitString(line);
				double x = atof(coordinates[0].c_str());
				double y = atof(coordinates[1].c_str());
				printf("Currently hardcoded in for traffic_lights\n");
				double z = atof(coordinates[2].c_str()) + 0.5;

				true_locations.push_back(R3Point(x, y, z));
			}

			truth.close();
		}
		printf("Arbitrary height and width, only works for traffic_lights\n");
		double object_width = 0.35;
		double object_height = 1;
		// go through every image
		for (int ir = 0; ir < scene->NRuns(); ir++) {
			GSVRun *run = scene->Run(ir);
			for (int is = 0, panorama_number = 0; is < scene->NSegments(); is++) {
				GSVSegment *segment = run->Segment(is);
				for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
					GSVPanorama *panorama = segment->Panorama(ip);
					for (int ii = 0; ii < panorama->NImages(); ii++) {
						GSVImage *image = panorama->Image(ii);

						// output filename
						char outputFilename[4096];
						sprintf(outputFilename, "proj_truth/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

						ofstream output;
						output.open(outputFilename, std::ofstream::trunc);

						if (!output.is_open()) { printf("Unable to write to %s\n", outputFilename); continue; }

						// project all of the true locations into this image
						for (unsigned int it = 0; it < true_locations.size(); it++) {
							R2Point position = image->UndistortedPosition(true_locations[it]);

							// make sure position falls within the image
							if (position.X() < 0 || position.X() > GSV_IMAGE_WIDTH) continue;
							if (position.Y() < 0 || position.Y() > GSV_IMAGE_HEIGHT) continue;

							// see how large the bounding box should be
							GSVCamera *camera = image->Camera();
							RNAngle xfov = 0.5 * camera->XFov();
							RNAngle yfov = 0.5 * camera->YFov();
							GSVPose pose = image->Pose();
							
							R3Point viewpoint = pose.Viewpoint();
							R3Vector towards = pose.Towards();
							R3Vector point = true_locations[it] - viewpoint;

							// project point onto towards for distance
							double distance = towards.Dot(point);
							
							double pixel_width = object_width / (2 * distance * tan(xfov)) * GSV_IMAGE_WIDTH;
							double pixel_height = object_height / (2 * distance * tan(yfov)) * GSV_IMAGE_HEIGHT;

							if (pixel_width < 20 || pixel_height < 20) continue;

							int width_radius = pixel_width / 2;
							int height_radius = pixel_height / 2;

							int xmin = position.X() - width_radius;
							int xmax = position.X() + width_radius;
							int ymin = GSV_IMAGE_HEIGHT - (position.Y() + height_radius);
							int ymax = GSV_IMAGE_HEIGHT - (position.Y() - height_radius);

							output << xmin;
							output << ',';
							output << ymin;
							output << ',';
							output << xmax;
							output << ',';
							output << ymax;
							output << ',';
							output << it;
							output << '\n';
						}
						output.close();
					}
				}
			}
		}
	}

	// Return success 
	return 0;
}