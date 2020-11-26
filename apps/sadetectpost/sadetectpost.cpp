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

static const char *detection_image_names[] = {
	"Detection", "Probabilities", "Visibility"
};

static int num_detection_image_names = sizeof(detection_image_names) / sizeof(const char *);

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
		fprintf(stderr, "Usage: sadetectpost input_scene [-v] [-split] <split>\n");
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
				RNScalar detection_zoom_factor = 1;

				// open up the position X grid
				char image_name[4096];
				sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_SA_Scanline.grd", run->Name(), is, ia);
				R2Grid scanline_grid;
				scanline_grid.Read(image_name);

				int xres = (int) (detection_zoom_factor * scanline_grid.XResolution());
				int yres = (int) (detection_zoom_factor * scanline_grid.YResolution());

				R2Affine world_to_grid = R2identity_affine;
				world_to_grid.Scale(detection_zoom_factor);
				world_to_grid.Transform(scanline_grid.WorldToGridTransformation());

				// create new grids
				for (int i = 0; i < num_detection_image_names; i++) {
					const char *grid_type = detection_image_names[i];

					R2Grid grid(xres, yres, world_to_grid);
					grid.Clear(R2_GRID_UNKNOWN_VALUE);

					// get spacing
					int spacing = scanline_grid.XResolution() / split;

					for (int x = 0; x < xres; x += spacing) {
						int xmin = x;
						int xmax = x + spacing;
						if (xmax > xres) xmax = xres;

						R2Grid smaller_grid;
						sprintf(image_name, "gsv_data/laser_images/%s/detections/%02d_%02d_SA_%06d_%06d_%s.grd", run->Name(), is, ia, xmin, xmax, grid_type);
						smaller_grid.Read(image_name);

						for (int xp = 0; xp < xmax - xmin; xp++) {
							for (int yp = 0; yp < yres; yp++) {
								grid.SetGridValue(x + xp, yp, smaller_grid.GridValue(xp, yp));
							}
						}
					}

					// write this grid 
					sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_SA_%s.grd", run->Name(), is, ia, grid_type);
					grid.Write(image_name);
				}
			}
		}
	}
	setup.close();

	// Return success 
	return 0;
}

















