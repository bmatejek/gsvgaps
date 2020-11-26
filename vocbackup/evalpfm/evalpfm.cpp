////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments
static int print_verbose = 0;
static const char *object = NULL;
static char *input_scene_name = NULL;
static char *unit = "SA";
static double max_allowable_distance = 1.0;

// detection object class
class Detection {
public:
	Detection(R3Point p, double score, int closest_match, double closest_distance);

public:
	R3Point p;
	double score;
	int closest_match;
	double closest_distance;
};

Detection::Detection(R3Point p, double score, int closest_match, double closest_distance) :
p(p),
score(score),
closest_match(closest_match),
closest_distance(closest_distance)
{
}

bool DetectionSort(Detection *a, Detection *b) {
	return a->score > b->score;
}

// Global variables
static GSVScene *scene = NULL;
static double epsilon = 10e-4;
static vector<Detection *> predictions = vector<Detection *>();
static vector<R3Point> true_locations = vector<R3Point>();

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *ReadScene(const char *filename)
{
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

////////////////////////////////////////////////////////////////////////
// Grid Reading Functions
////////////////////////////////////////////////////////////////////////

static int ReadGrids(void) {
	RNTime start_time;
	start_time.Read();
	if (print_verbose) {
		printf("Reading grids ...\n");
		fflush(stdout);
	}

	// create a list of prediction points
	int number_scans = 3, zero_probability = 0;
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); is++) {
			for (int i = 0; i < number_scans; i++) {
				// never use the middle scan
				if (i == 1) continue;
				char image_name[4096];
				R2Grid x_position, y_position, z_position, detection;
				
				// read in all of the laser images
				sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_PositionX.grd", run->Name(), is, i, unit);
				if (!x_position.Read(image_name)) return 0;
				sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_PositionY.grd", run->Name(), is, i, unit);
				if (!y_position.Read(image_name)) return 0;
				sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_PositionZ.grd", run->Name(), is, i, unit);
				if (!z_position.Read(image_name)) return 0;
				sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_%s_Mean.grd", run->Name(), is, i, unit, object);
				if (!detection.Read(image_name)) return 0;

				for (int px = 0; px < detection.XResolution(); px++) {
					for (int py = 0; py < detection.YResolution(); py++) {
						R3Point p = R3Point(x_position.GridValue(px, py), y_position.GridValue(px, py), z_position.GridValue(px, py));
						double score = detection.GridValue(px, py);
						if (score > epsilon) {
							int closest_match = 0;
							double closest_distance = DBL_MAX;

							// find the closest location
							for (unsigned int t = 0; t < true_locations.size(); t++) {
								if (R3Distance(p, true_locations[t]) < closest_distance) {
									closest_distance = R3Distance(p, true_locations[t]);
									closest_match = t;
								}
							}

							Detection *prediction = new Detection(p, score, closest_match, closest_distance);
							predictions.push_back(prediction);
						}
						else {
							zero_probability++;
						}
					}
				}
			}
		}
	}
	
	// Print statistics
	if (print_verbose) {
		printf("  Done in %.2f seconds\n", start_time.Elapsed());
		printf("  Number of predictions: %lu\n", predictions.size());
		printf("  Number of discarded locations: %d\n", zero_probability);
		fflush(stdout);
	}	
	
	start_time.Read();
	if (print_verbose) {
		printf("Evaluating detection algorithm ...\n");
		fflush(stdout);
	}

	sort(predictions.begin(), predictions.end(), DetectionSort);
	int false_positives = 0;
	int matched = 0;
	bool *match = new bool[true_locations.size()];
	for (int i = 0; i < true_locations.size(); i++) {
		match[i] = false;
	}
	for (int t = 0; t < predictions.size(); t++) {
		if (predictions[t]->closest_distance > max_allowable_distance) {
			false_positives++;
		}
		else {
			if (!match[predictions[t]->closest_match]) {
				match[predictions[t]->closest_match] = true;
				matched++;
				printf("%f\n", matched / (double)(matched + false_positives));
			}
		}
	}

	printf("Recall: %f\n", matched / (double)true_locations.size());

	// return success
	return 1;
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
	fprintf(stderr, "Usage: evalpfm [input_scene] -object <object> [-max_distance <max_allowable_distance>]\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-object")) { argv++; argc--; object = *argv; }
			else if (!strcmp(*argv, "-DA")) { unit = "DA"; }
			else if (!strcmp(*argv, "-max_distance")) { argv++; argc--; max_allowable_distance = atof(*argv); }
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) { input_scene_name = *argv; }
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// Check all input
	if (!input_scene_name) return Usage();
	if (!object) return Usage();

	// Return OK status 
	return 1;
}


static int GetTruthPoints(GSVScene *scene) {
	// get all truth data points
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);

		char truthFilename[128];
		sprintf(truthFilename, "truth/%s/%s/truth.txt", object, run->Name());
		ifstream truth;
		truth.open(truthFilename);

		if (!truth.is_open()) { fprintf(stderr, "Failed to read truth filename: %s\n", truthFilename); return 1; }

		string line;
		while (getline(truth, line)) {
			vector<string> location = splitString(line);
			double x = atof(location[0].c_str());
			double y = atof(location[1].c_str());
			printf("Currently hard coded in for stop signs\n");
			double z = atof(location[2].c_str()) + 0.375;

			true_locations.push_back(R3Point(x, y, z));
		}

		truth.close();
	}
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// read in the scene
	scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);

	// read all of the truth points
	if (!GetTruthPoints(scene)) exit(-1);

	// read in all of the grids
	if (!ReadGrids()) exit(-1);

	// Return success 
	return 0;
}