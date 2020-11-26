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

// Program arguments
static const char *input_scene_name = NULL;
static const char *output_directory = NULL;
static int print_verbose = 0;
static int split = 200;							// number of ways to split total

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

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-split")) { argv++; argc--; split = atoi(*argv); }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
			argv++; argc--;
		}
	}

	// Check scene name
	if (!input_scene_name) {
		fprintf(stderr, "Usage: sadetectsetup input_scene [-v] [-split] <split>\n");
		return FALSE;
	}

	// Check output image directory
	if (!output_directory) {
		output_directory = "scripts";
	}

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

	// open up the temporary output file
	char filename[4096];
	sprintf(filename, "%s/setup.txt", output_directory);
	ofstream setup;
	setup.open(filename);
	if (!setup.is_open()) { fprintf(stderr, "Error: failed to open %s\n", filename); return 0; }

	// go through each segment and scan
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ia = 0; ia < segment->NScans(); ia++) {
				if (ia == 1) continue;
				// open up SA_Scanline
				R2Grid scanline_grid;
				sprintf(filename, "gsv_data/laser_images/%s/%02d_%02d_SA_Scanline.grd", run->Name(), is, ia);
				scanline_grid.Read(filename);
				setup << ir;
				setup << ',';
				setup << is;
				setup << ',';
				setup << ia;
				setup << ',';
				setup << scanline_grid.XResolution();
				setup << ',';
				setup << scanline_grid.XResolution() / split;
				setup << '\n';
			}
		}
	}
	setup.close();

	// Return success 
	return 0;
}

















