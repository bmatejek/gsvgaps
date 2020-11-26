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
static const char *output_image_directory = NULL;
static int print_verbose = 0;
static int print_debug = 0;

static char *object = NULL;
static double high_threshold = 1.0;
static double low_threshold = -0.5;
static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;

// XY, DA, DH parameters

static RNLength SA_viewpoint_spacing = 0.01;
static RNLength DA_image_spacing = 0.1;

static double min_frustum_depth = 1;
static double max_frustum_depth = 1000;
static double max_viewpoint_distance = -1;



////////////////////////////////////////////////////////////////////////
// Grid names
////////////////////////////////////////////////////////////////////////

const char *object_names[] = {
	"Mean", "Visibility"
};

int num_object_names = sizeof(object_names) / sizeof(const char *);

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

static double ConvertScore(double score) {
	if (score < low_threshold) return 0.0;
	else if (score > high_threshold) return 1.0;
	else return (score - low_threshold) / (high_threshold - low_threshold);
}

////////////////////////////////////////////////////////////////////////
// Utility types and functions
////////////////////////////////////////////////////////////////////////

static void InterpolateMissingColumns(R2Grid& grid)
{
	// Initialize previous column with values
	int ix0 = -1;

	// Consider every column
	for (int ix1 = 0; ix1 < grid.XResolution(); ix1++) {
		// Check if column has values
		RNBoolean found_value = FALSE;
		for (int iy = 0; iy < grid.YResolution(); iy++) {
			if (grid.GridValue(ix1, iy) != R2_GRID_UNKNOWN_VALUE) {
				found_value = TRUE;
				break;
			}
		}

		// Skip column, if has no values
		if (!found_value) continue;

		// Iterpolate values in skipped columns
		if ((ix0 >= 0) && (ix0 < ix1 - 1)) {
			for (int iy = 0; iy < grid.YResolution(); iy++) {
				RNScalar value0 = grid.GridValue(ix0, iy);
				if (value0 == R2_GRID_UNKNOWN_VALUE) continue;
				RNScalar value1 = grid.GridValue(ix1, iy);
				if (value1 == R2_GRID_UNKNOWN_VALUE) continue;
				for (int ix = ix0 + 1; ix < ix1; ix++) {
					RNScalar t = (double)(ix - ix0) / (double)(ix1 - ix0);
					RNScalar value = (1 - t)*value0 + t*value1;
					grid.SetGridValue(ix, iy, value);
				}
			}
		}

		// Remember last column with values
		ix0 = ix1;
	}
}

static int DetectObjects(GSVScan *scan, const char *output_image_directory, const char *image_parameterization) {
	// Parameters
	RNScalar zoom_factor = 1;

	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	
	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;
	GSVScene *scene = run->Scene();
	const char *cache_directory = scene->CacheDataDirectoryName();
	char image_name[4096];

	// padding is the number of panoramas before this segment
	int padding = 0;
	for (int is = 0; is < segment_index; is++) {
		padding += run->Segment(is)->NPanoramas();
	}

	// Print message
	if (print_debug) {
		printf("    Creating object detection images ...\n");
		fflush(stdout);
	}

	// Read scanline image
	R2Grid scanline_grid, timestamp_grid;
	sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_Scanline.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!scanline_grid.Read(image_name)) return 0;

	// Read position images
	R2Grid position_x_grid, position_y_grid, position_z_grid;
	sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionX.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!position_x_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionY.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!position_y_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionZ.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!position_z_grid.Read(image_name)) return 0;

	// Compute world to grid transformation
	int xres = (int)(zoom_factor * scanline_grid.XResolution());
	int yres = (int)(zoom_factor * scanline_grid.YResolution());

	// Create object images for scan
	R2Affine world_to_grid = R2identity_affine;
	world_to_grid.Scale(zoom_factor);
	world_to_grid.Transform(scanline_grid.WorldToGridTransformation());
	R2Grid average(xres, yres, world_to_grid);
	average.Clear(0);
	R2Grid visibility(average);

	// initialize non scans  parts to R2_GRID_UNKNOWN_VALUE
	for (int ix = 0; ix < xres; ix++) {
		for (int iy = 0; iy < yres; iy++) {
			// find this pixel
			int px = (int)(ix / zoom_factor + 0.5);
			if (px >= position_x_grid.XResolution()) continue;
			int py = (int)(iy / zoom_factor + 0.5);
			if (py >= position_x_grid.YResolution()) continue;

			visibility.SetGridValue(px, py, R2_GRID_UNKNOWN_VALUE);

			// get (x, y, z) world coordinates for this pixel
			RNScalar x = position_x_grid.GridValue(px, py);
			if (x == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar y = position_y_grid.GridValue(px, py);
			if (y == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar z = position_z_grid.GridValue(px, py);
			if (z == R2_GRID_UNKNOWN_VALUE) continue;
			visibility.SetGridValue(px, py, 0.0);
		}
	}

	int ntapestries[3] = { 8, 0, 8 };
	int tapestries[3][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 0 }, { 0, 1, 2, 3, 4, 5, 6, 7 } };
	for (int k = 0; k < ntapestries[scan_index]; k++) {
		int it = tapestries[scan_index][k];
		GSVTapestry *tapestry = segment->Tapestry(it);
		GSVCamera *camera = tapestry->Camera();
		if (!camera) continue;

		// Consider every image in tapestry
		for (int ii = 0; ii < tapestry->NImages(); ii++) {
			GSVImage *image = tapestry->Image(ii);
			// Get image info
			const GSVPose& image_pose = image->Pose();
			RNScalar image_timestamp = image->Timestamp();
			const R3Point& image_viewpoint = image_pose.Viewpoint();
			R3Frustum image_frustum(image_viewpoint, image_pose.Towards(), image_pose.Up(),	camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);

			// get the bounding box predictions
			char probabilityFilename[128];
			sprintf(probabilityFilename, "voc_predictions/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), segment_index, padding + image->TapestryIndex(), image->PanoramaIndex());
			ifstream probabilities;
			probabilities.open(probabilityFilename);
			if (!probabilities.is_open()) { fprintf(stderr, "Failed to open %s\n", probabilityFilename); return 0; }

			if (print_verbose) printf("Reading %s\n", probabilityFilename);

			// read the predictions into a grid
			R2Grid *image_probabilities = new R2Grid(GSV_IMAGE_WIDTH, GSV_IMAGE_HEIGHT);
			string line;
			while (getline(probabilities, line)) {
				vector<string> prediction = splitString(line);

				// get coordinates of box
				double xmin = atof(prediction[0].c_str());
				double ymax = image->Height() - atof(prediction[1].c_str());
				double xmax = atof(prediction[2].c_str());
				double ymin = image->Height() - atof(prediction[3].c_str());
				double score = atof(prediction[5].c_str());

				double probability = ConvertScore(score);

				// add in this probability 
				for (int x = floor(xmin); x <= ceil(xmax); x++) {
					for (int y = floor(ymin); y <= ceil(ymax); y++) {
						double old_probability = image_probabilities->GridValue(x, y);
						if (old_probability < probability) {
							image_probabilities->SetGridValue(x, y, probability);
						}
					}
				}

			}
			probabilities.close();

			// Consider every grid pixel
			for (int ix = 0; ix < xres; ix++) {
				for (int iy = 0; iy < yres; iy++) {
					// find this pixel
					int px = (int)(ix / zoom_factor + 0.5);
					if (px >= position_x_grid.XResolution()) continue;
					int py = (int)(iy / zoom_factor + 0.5);
					if (py >= position_x_grid.YResolution()) continue;
					
					// get (x, y, z) world coordinates for this pixel
					RNScalar x = position_x_grid.GridValue(px, py);
					if (x == R2_GRID_UNKNOWN_VALUE) continue;
					RNScalar y = position_y_grid.GridValue(px, py);
					if (y == R2_GRID_UNKNOWN_VALUE) continue;
					RNScalar z = position_z_grid.GridValue(px, py);
					if (z == R2_GRID_UNKNOWN_VALUE) continue;

					// make sure there is a scanline here
					RNScalar scanline_value = scanline_grid.GridValue(px, py);
					//if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
					int scanline_index = (int)(scanline_value + 0.5);

					// get scanline
					GSVScanline *scanline = scan->Scanline(scanline_index);
					R3Point position(x, y, z);

					// Check if scanline viewpoint is within range
					const GSVPose& scanline_pose = scanline->Pose();
					if (max_viewpoint_distance > 0) {
						if (R3Distance(scanline_pose.Viewpoint(), image_viewpoint) > max_viewpoint_distance) continue;
					}

					// Check if point is within image view frustum
					if (!image_frustum.Intersects(position)) continue;

					// Retrieve estimate of point from undistorted image
					R2Point image_position = image->UndistortedPosition(position);

					// make sure this is a valid image location
					if (image_position.X() < 0) continue;
					if (image_position.Y() < 0) continue;
					if (image_position.X() >= image_probabilities->XResolution()) continue;
					if (image_position.Y() >= image_probabilities->YResolution()) { continue; }
					int imageX = (int)round(image_position.X());
					int imageY = (int)round(image_position.Y());

					double probability = image_probabilities->GridValue(imageX, imageY);

					// Update scan grid
					if (probability > 0) {
						average.SetGridValue(ix, iy, average.GridValue(ix, iy) + probability);
					}
					visibility.SetGridValue(ix, iy, visibility.GridValue(ix, iy) + 1);
				}
			}
			delete image_probabilities;
		}
	}

	// update average
	double epsilon = 10e-1;
	for (int px = 0; px < average.XResolution(); px++) {
		for (int py = 0; py < average.YResolution(); py++) {
			if (visibility.GridValue(px, py) > epsilon) {
				average.SetGridValue(px, py, average.GridValue(px, py) / visibility.GridValue(px, py));
			}
			else {
				average.SetGridValue(px, py, R2_GRID_UNKNOWN_VALUE);
			}
		}
	}

	// Write scan grids 
	sprintf(image_name, "%s/%s/%02d_%02d_%s_%s_Mean.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, object);
	average.Write(image_name);

	sprintf(image_name, "%s/%s/%02d_%02d_%s_Visibility.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	visibility.Write(image_name);

	// Return success
	return 1;
}


////////////////////////////////////////////////////////////////////////
// SA images (Scanline vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteSAImages(GSVScan *scan, const char *output_image_directory) {
	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();

	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Print message
	if (print_verbose) {
		printf("  Creating SA images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
		fflush(stdout);
	}

	// Capture SA images
	if (!DetectObjects(scan, output_image_directory, "SA")) return 0;

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// DA images (Distance vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteDAImage(GSVScan *scan, const char *output_image_directory, const R2Grid& sa_travel_distance_grid, const R2Grid& sa_viewpoint_distance_grid, const char *grid_name) {
	if (!strcmp(grid_name, "Mean")) {
		grid_name = new char[128];
		grid_name = (string(object) + string("_") + string("Mean")).c_str();
	}

	// Get convenient variables
	if (!scan) return 0;
	if (scan->NScanlines() == 0) return 0;
	int scan_index = scan->SegmentIndex();

	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Skip if already done
	char da_name[4096];
	sprintf(da_name, "%s/%s/%02d_%02d_DA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
	//if (RNFileExists(da_name)) return 1;

	// Read SA image
	R2Grid sa_grid;
	char sa_name[4096];
	sprintf(sa_name, "%s/%s/%02d_%02d_SA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
	sa_grid.Read(sa_name);

	// Compute grid resolution
	RNScalar total_travel_distance = sa_travel_distance_grid.Maximum();
	if (total_travel_distance == 0) return 0;
	int xres = (int)(total_travel_distance / DA_image_spacing + 0.5);
	if (xres == 0) return 0;
	int yres = 180;

	// Compute world to grid transformation
	R2Affine world_to_grid = R2identity_affine;
	world_to_grid.XScale(1.0 / DA_image_spacing);

	// Create DA image
	R2Grid da_grid(xres, yres, world_to_grid);
	da_grid.Clear(R2_GRID_UNKNOWN_VALUE);

	// Create DA distance image
	R2Grid da_viewpoint_distance_grid(xres, yres, world_to_grid);
	da_viewpoint_distance_grid.Clear(FLT_MAX);

	// Fill DA image
	for (int i = 0; i < sa_grid.XResolution(); i++) {
		for (int iy = 0; iy < sa_grid.YResolution(); iy++) {
			RNScalar value = sa_grid.GridValue(i, iy);
			if (value == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar viewpoint_distance = sa_viewpoint_distance_grid.GridValue(i, iy);
			if (viewpoint_distance == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar travel_distance = sa_travel_distance_grid.GridValue(i, iy);
			if (travel_distance == R2_GRID_UNKNOWN_VALUE) continue;
			int ix = (int)(xres * travel_distance / total_travel_distance + 0.5);
			if ((ix < 0) || (ix >= xres)) continue;
			RNScalar old_viewpoint_distance = da_viewpoint_distance_grid(ix, iy);
			if (viewpoint_distance < old_viewpoint_distance)  {
				da_viewpoint_distance_grid.SetGridValue(ix, iy, viewpoint_distance);
				da_grid.SetGridValue(ix, iy, value);
			}
		}
	}

	// Fill holes in DA image
	InterpolateMissingColumns(da_grid);

	// Write DA image
	da_grid.Write(da_name);

	// Return success
	return 1;
}

static int WriteDAImages(GSVScan *scan, const char *output_image_directory) {
	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	
	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Print message
	if (print_verbose) {
		printf("  Creating DA images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
		fflush(stdout);
	}

	// Read SA scanline and distance images
	char sa_travel_distance_name[4096];
	char sa_viewpoint_distance_name[4096];
	R2Grid sa_travel_distance_grid;
	R2Grid sa_viewpoint_distance_grid;
	sprintf(sa_travel_distance_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	sprintf(sa_viewpoint_distance_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	if (!RNFileExists(sa_travel_distance_name) && (!WriteSAImages(scan, output_image_directory))) return 0;
	sa_travel_distance_grid.Read(sa_travel_distance_name);
	sa_viewpoint_distance_grid.Read(sa_viewpoint_distance_name);

	// Write DA images
	for (int i = 0; i < num_object_names; i++) {
		if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, object_names[i])) return 0;
	}

	// Return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Top-level image computation function
////////////////////////////////////////////////////////////////////////

static int WriteImages(GSVScene *scene, const char *output_image_directory) {
	// Start statistics
	RNTime start_time;
	start_time.Read();
	if (print_verbose) {
		printf("Creating images ...\n");
		fflush(stdout);
	}

	// Write scan images
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ia = 0; ia < segment->NScans(); ia++) {
				GSVScan *scan = segment->Scan(ia);
				if (!WriteSAImages(scan, output_image_directory)) return 0;
				if (!WriteDAImages(scan, output_image_directory)) return 0;
			}
		}
	}

	// Print statistics
	if (print_verbose) {
		printf("  Done in %.2f seconds\n", start_time.Elapsed());
		fflush(stdout);
	}

	// Return success
	return 1;
}




////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage() {
	fprintf(stderr, "Usage: dadetect input_scene [output_directory] -object <object> [options]\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-debug")) { print_debug = 1; }
			else if (!strcmp(*argv, "-object")) { argv++; argc--; object = *argv; }
			else if (!strcmp(*argv, "-max_distance")) { argv++; argc--; max_viewpoint_distance = atof(*argv); }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else if (!output_image_directory) output_image_directory = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// Check all input
	if (!input_scene_name) return Usage();
	if (!object) return Usage();

	// Check output image directory
	if (!output_image_directory) {
		output_image_directory = "gsv_data/laser_images";
	}

	// Check verbosity 
	if (print_debug) {
		print_verbose = 1;
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

	// Create output image directories
	char cmd[1024];
	sprintf(cmd, "mkdir -p %s", output_image_directory);
	system(cmd);
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		sprintf(cmd, "mkdir -p %s/%s", output_image_directory, run->Name());
		system(cmd);
	}

	// Write images
	if (!WriteImages(scene, output_image_directory)) exit(-1);

	// Return success 
	return 0;
}