////////////////////////////////////////////////////////////////////////
// Include Files
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
// Program Arguments
////////////////////////////////////////////////////////////////////////

static int print_verbose = 0;
static GSVScene *scene = NULL;
static char *input_scene_name = NULL;
static char *object = NULL;

class ImagePrediction {
public:
	ImagePrediction(R2Box prediction, double score, bool match, int matched_with);

public:
	R2Box prediction;
	double score;
	bool match;
	int matched_with;
};

ImagePrediction::ImagePrediction(R2Box prediction, double score, bool match, int matched_with) :
prediction(prediction),
score(score),
match(match),
matched_with(matched_with)
{
}

class TruePrediction {
public:
	TruePrediction(R2Box truth, int index);

public:
	R2Box truth;
	int index;
};

TruePrediction::TruePrediction(R2Box truth, int index) :
truth(truth),
index(index)
{
}

static vector<ImagePrediction *> predictions = vector<ImagePrediction *>();
static vector<R3Point> true_locations = vector<R3Point>();

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
vector<string> splitString(string line, char delim = ',') {
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
			double z = atof(location[2].c_str());

			true_locations.push_back(R3Point(x, y, z));
		}

		truth.close();
	}
	return 1;
}


bool PredictionSort(ImagePrediction *a, ImagePrediction *b) {
	return a->score > b->score;
}


static int ReadAll(void) {
	// open up predictions file
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0, panorama_number = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ii++) {
					// read in the projected truth for this image
					vector<TruePrediction *> proj_boxes = vector<TruePrediction *>();

					char projFilename[4096];
					sprintf(projFilename, "proj_truth/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

					ifstream truth;
					truth.open(projFilename);
					if (!truth.is_open()) { fprintf(stderr, "Failed to read %s\n", projFilename); return 0; }

					string line;
					while (getline(truth, line)) {
						vector<string> coordinates = splitString(line);
						if (coordinates.size() != 5) { fprintf(stderr, "Failed to read %s\n", projFilename); return 0; }
						double xmin = atof(coordinates[0].c_str());
						double ymin = atof(coordinates[1].c_str());
						double xmax = atof(coordinates[2].c_str());
						double ymax = atof(coordinates[3].c_str());

						int index = atoi(coordinates[4].c_str());

						R2Box box = R2Box(xmin, ymin, xmax, ymax);
						proj_boxes.push_back(new TruePrediction(box, index));
					}

					char predFilename[4096];
					sprintf(predFilename, "voc_predictions/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

					ifstream pred;
					pred.open(predFilename);
					if (!pred.is_open()) { fprintf(stderr, "Failed to open %s\n", predFilename); return 0; }

					while (getline(pred, line)) {
						vector<string> coordinates = splitString(line);
						if (coordinates.size() != 6) { fprintf(stderr, "Failed to read %s\n", predFilename); return 0; }

						double xmin = atof(coordinates[0].c_str());
						double ymin = atof(coordinates[1].c_str());
						double xmax = atof(coordinates[2].c_str());
						double ymax = atof(coordinates[3].c_str());

						R2Box box = R2Box(xmin, ymin, xmax, ymax);
						double score = atof(coordinates[5].c_str());
						bool match = false;
						int matched_with = -1;

						// go through every true bounding box in this image to see if there is a match
						for (unsigned int it = 0; it < proj_boxes.size(); it++) {
							R2Box union_truth = proj_boxes[it]->truth;
							R2Box intersection_truth = proj_boxes[it]->truth;

							union_truth.Union(box);
							intersection_truth.Intersect(box);

							// see if the intersection over union is greater than 0.5
							if (intersection_truth.Area() / union_truth.Area() > 0.10) {
								match = true;
								matched_with = proj_boxes[it]->index;
								break;
							}
						}

						ImagePrediction *image_prediction = new ImagePrediction(box, score, match, matched_with);
						predictions.push_back(image_prediction);
					}
					pred.close();
				}
			}
		}
	}

	// sort the  image predictions
	sort(predictions.begin(), predictions.end(), PredictionSort);
	int false_positives = 0;
	int matched = 0;
	bool *match = new bool[true_locations.size()];
	for (unsigned int i = 0; i < true_locations.size(); i++) {
		match[i] = false;
	}
	for (unsigned int t = 0; t < predictions.size(); t++) {
		if (!predictions[t]->match) {
			false_positives++;
		}
		else {
			if (!match[predictions[t]->matched_with]) {
				match[predictions[t]->matched_with] = true;
				matched++;
				printf("%f\n", matched / (double)(matched + false_positives));
			}
		}
	}
	
	printf("Recall: %f\n", matched / (double)true_locations.size());

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage() {
	fprintf(stderr, "Usage: evalpred [input_scene] -object <object>\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-object")) { argv++; argc--; object = *argv; }
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

	// read in all predictions and truth
	if (!ReadAll()) exit(-1);

	// return OK status 
	return 1;
}
