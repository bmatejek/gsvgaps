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
static const char *output_image_directory = NULL;
static int print_verbose = 0;
static int xscan_min = -1;
static int xscan_max = -1;
static int input_run = -1;
static int input_segment = -1;
static int input_scan = -1;

// Depth and color parameters
static int low_scans_occluded_by_car = 25;
static double min_frustum_depth = 1;
static double max_frustum_depth = 1000;
static double max_viewpoint_distance = 100;

// Global variables
static int GSV_IMAGE_WIDTH = 1936;
static int GSV_IMAGE_HEIGHT = 2592;
static double low_threshold = 0;
static double high_threshold = 3.0;

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
// Object Detection Images
////////////////////////////////////////////////////////////////////////

static int WriteDetectionImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization) {
	// open up the objects file
	char objectListFilename[4096];
	sprintf(objectListFilename, "rcnn/objects.txt");

	ifstream objects;
	objects.open(objectListFilename);
	if (!objects.is_open()) { fprintf(stderr, "Failed to read objects file at %s\n", objectListFilename); return 0; }

	map<string, int> object_mapping;
	string line;
	while (getline(objects, line)) {
		vector<string> object = splitString(line, ' ');
		object_mapping.insert(pair<string, int>(object[1], atoi(object[0].c_str())));
	}
	objects.close();
	// Parameters
	RNScalar detection_zoom_factor = 1;

	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
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
	if (print_verbose) {
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
	if (xscan_min == -1) xscan_min = 0;
	if (xscan_max == -1) xscan_max = scanline_grid.XResolution();
	int xres = (int)(detection_zoom_factor * (xscan_max - xscan_min));
	int yres = (int)(detection_zoom_factor * scanline_grid.YResolution());

	// Create detection images for scan
	R2Affine world_to_grid = R2identity_affine;
	world_to_grid.Scale(detection_zoom_factor);
	world_to_grid.Transform(scanline_grid.WorldToGridTransformation());
	R2Grid detection_grid(xres, yres, world_to_grid);
	detection_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid visibility_grid(detection_grid);
	R2Grid probability_grid(detection_grid);

	// initialize non scan parts to R2_GRID_UNKNOWN_VALUE
	for (int ix = xscan_min, ixp = 0; ix < xscan_max; ++ix, ++ixp) {
		for (int iy = 0; iy < yres; iy++) {
			// find this pixel
			int px = (int)(ix / detection_zoom_factor + 0.5);
			if (px >= position_x_grid.XResolution()) continue;
			int py = (int)(iy / detection_zoom_factor + 0.5);
			if (py >= position_x_grid.YResolution()) continue;

			// default value for visibility grid is R2_GRID_UNKNOWN_VALUE if there is no corresponding (x, y, z) location
			visibility_grid.SetGridValue(ixp, py, R2_GRID_UNKNOWN_VALUE);

			// otherwise the visibility grid gets the value 0 to start
			RNScalar x = position_x_grid.GridValue(px, py);
			if (x == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar y = position_y_grid.GridValue(px, py);
			if (y == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar z = position_z_grid.GridValue(px, py);
			if (z == R2_GRID_UNKNOWN_VALUE) continue;
			if (iy < low_scans_occluded_by_car) visibility_grid.SetGridValue(ixp, py, R2_GRID_UNKNOWN_VALUE);
			else visibility_grid.SetGridValue(ixp, py, 0.0);

		}
	}
	// Consider tapestries relevant to scan
	int ntapestries[3] = { 8, 0, 8 };
	int tapestries[3][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 }, { 0 }, { 0, 1, 2, 3, 4, 5, 6, 7 } };

	// keep an array of open files
	ifstream **rcnn_predictions = new ifstream *[segment->NImages()];
	if (scan_index == 0) return 1;
	for (int k = 0, c = 0; k < ntapestries[scan_index]; k++) {
		int it = tapestries[scan_index][k];
		GSVTapestry *tapestry = segment->Tapestry(it);
		for (int ii = 0; ii < tapestry->NImages(); ++ii, ++c) {
			GSVImage *image = tapestry->Image(ii);
			char rcnnFilename[4096];
			sprintf(rcnnFilename, "rcnn/%s/%02d_%06d_%02d_UndistortedImage.txt", run->Name(), segment_index, padding + image->TapestryIndex(), image->PanoramaIndex());

			// open up this set of predictions
			ifstream *prediction = new ifstream();
			prediction->open(rcnnFilename);

			// send error if this prediction does not open
			if (!prediction->is_open()) { fprintf(stderr, "Failed to open %s\n", rcnnFilename); return 0; }
			rcnn_predictions[c] = prediction;
		}
	}

	// consider every pixel in image
	for (int ix = xscan_min, ixp = 0; ix < xscan_max; ++ix, ++ixp) {
		for (int iy = 0; iy < yres; iy++) {
			// find this pixel
			int px = (int)(ix / detection_zoom_factor + 0.5);
			if (px >= position_x_grid.XResolution()) continue;
			int py = (int)(iy / detection_zoom_factor + 0.5);
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
			if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
			int scanline_index = (int)(scanline_value + 0.5);

			// get scanline
			GSVScanline *scanline = scan->Scanline(scanline_index);
			R3Point position(x, y, z);
			// consider every object
			double *probabilities = new double[object_mapping.size()];
			// initialize to 0
			for (int i = 0; i < (int)object_mapping.size(); i++) {
				probabilities[i] = 0;
			}

			for (int k = 0, c = 0; k < ntapestries[scan_index]; k++) {
				int it = tapestries[scan_index][k];
				GSVTapestry *tapestry = segment->Tapestry(it);
				GSVCamera *camera = tapestry->Camera();
				if (!camera) continue;

				// Consider every image in tapestry
				for (int ii = 0; ii < tapestry->NImages(); ++ii, ++c) {
					GSVImage *image = tapestry->Image(ii);

					// Get image info
					const GSVPose& image_pose = image->Pose();
					const R3Point& image_viewpoint = image_pose.Viewpoint();

					// Check if scanline viewpoint is within range
					const GSVPose& scanline_pose = scanline->Pose();

					// Check if point is within image view frustum
					R3Frustum image_frustum(image_viewpoint, image_pose.Towards(), image_pose.Up(), camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);
					if (!image_frustum.Intersects(position)) continue;

					// Retrieve estimate of point from undistorted image
					R2Point image_position = image->UndistortedPosition(position);

					// make sure this is a valid image location
					if (image_position.X() < 0) continue;
					if (image_position.Y() < 0) continue;
					if (image_position.X() >= GSV_IMAGE_WIDTH) continue;
					if (image_position.Y() >= GSV_IMAGE_HEIGHT) continue;
					// read the predictions
					ifstream *prediction = rcnn_predictions[c];

					string line;
					while (getline(*prediction, line)) {
						vector<string> pred = splitString(line);
						// get coordinates of box
						string object = pred[0];
						double score = atof(pred[1].c_str());
						double xmin = atof(pred[2].c_str());
						double ymax = image->Height() - atof(pred[3].c_str());
						double xmax = atof(pred[4].c_str());
						double ymin = image->Height() - atof(pred[5].c_str());
						R2Box box = R2Box(xmin, ymin, xmax, ymax);

						// see if this box contains the image point
						if (!(image_position.X() > xmin && image_position.X() < xmax && image_position.Y() > ymin && image_position.Y() < ymax)) continue;

						if (score > 0.0 && R3Distance(R3Point(x, y, z), R3Point(93.841930, 827.364207, -30.804808)) < 1.0) {
							printf("%s: %f\n", object.c_str(), score);
						}

						// only include this point if it is sufficiently far from the edge
						double midOffsetY = abs(image_position.Y() - (ymax - ymin) / 2 - ymin);
						if (midOffsetY / (ymax - ymin) > 0.1) continue;
						double midOffsetX = abs(image_position.X() - (xmax - xmin) / 2 - xmin);
						if (midOffsetX / (xmax - xmin) > 0.1) continue;

						// get key value
						int key = object_mapping.at(object);
						probabilities[key] += ConvertScore(score);
					}
					// reset the EOF flag
					prediction->clear();
					prediction->seekg(0, ios::beg);

					// set the values for this point
					visibility_grid.SetGridValue(ixp, iy, 1 + visibility_grid.GridValue(ixp, iy));
				}
			}
			if (visibility_grid.GridValue(ixp, iy) == 0) {
				detection_grid.SetGridValue(ixp, iy, R2_GRID_UNKNOWN_VALUE);
				probability_grid.SetGridValue(ixp, iy, R2_GRID_UNKNOWN_VALUE);
			}
			else {
				// find largest prediction
				int best_object = 0;
				double best_probability = DBL_MIN;
				double sum_probabilities = 0;
				for (int i = 0; i < (int)object_mapping.size(); i++) {
					sum_probabilities += probabilities[i];
					if (best_probability >= probabilities[i]) continue;
					best_probability = probabilities[i];
					best_object = i;
				}

				if (iy < low_scans_occluded_by_car) {
					probability_grid.SetGridValue(ixp, iy, R2_GRID_UNKNOWN_VALUE);
					detection_grid.SetGridValue(ixp, iy, R2_GRID_UNKNOWN_VALUE);
				}
				else if (best_probability == 0.0) {
					probability_grid.SetGridValue(ixp, iy, 0.0);
					detection_grid.SetGridValue(ixp, iy, R2_GRID_UNKNOWN_VALUE);
				}
				else {
					probability_grid.SetGridValue(ixp, iy, best_probability);
					detection_grid.SetGridValue(ixp, iy, best_object);
				}
			}
		}
	}

	// close rcnn files
	for (int k = 0, c = 0; k < ntapestries[scan_index]; k++) {
		int it = tapestries[scan_index][k];
		GSVTapestry *tapestry = segment->Tapestry(it);
		for (int ii = 0; ii < tapestry->NImages(); ++ii, ++c) {
			rcnn_predictions[c]->close();
		}
	}

	// write scan grids
	if (xscan_min == -1 && xscan_max == -1) {
		sprintf(image_name, "%s/%s/%02d_%02d_%s_%06d_%06d_Detection.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		detection_grid.Write(image_name);
		sprintf(image_name, "%s/%s/%02d_%02d_%s_%06d_%06d_Probabilities.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		probability_grid.Write(image_name);
		sprintf(image_name, "%s/%s/%02d_%02d_%s_%06d_%06d_Visibility.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		visibility_grid.Write(image_name);
	}
	else {
		sprintf(image_name, "%s/%s/detections/%02d_%02d_%s_%06d_%06d_Detection.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		detection_grid.Write(image_name);
		sprintf(image_name, "%s/%s/detections/%02d_%02d_%s_%06d_%06d_Probabilities.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		probability_grid.Write(image_name);
		sprintf(image_name, "%s/%s/detections/%02d_%02d_%s_%06d_%06d_Visibility.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, xscan_min, xscan_max);
		visibility_grid.Write(image_name);
	}

	// Return success
	return 1;
}

static int WriteSAImages(GSVScan *scan, const char *output_image_directory) {
	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
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
	if (!WriteDetectionImages(scan, output_image_directory, "SA")) return 0;

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
	if (input_run == -1 && input_segment == -1 && input_scan == -1) {
		for (int ir = 0; ir < scene->NRuns(); ir++) {
			GSVRun *run = scene->Run(ir);
			for (int is = 0; is < run->NSegments(); is++) {
				GSVSegment *segment = run->Segment(is);
				for (int ia = 0; ia < segment->NScans(); ia++) {
					// don't need middle scan
					if (ia == 1) continue;
					GSVScan *scan = segment->Scan(ia);
					if (!WriteSAImages(scan, output_image_directory)) return 0;
				}
			}
		}
	}
	else {
		GSVRun *run = scene->Run(input_run);
		GSVSegment *segment = run->Segment(input_segment);
		GSVScan *scan = segment->Scan(input_scan);
		if (!WriteSAImages(scan, output_image_directory)) return 0;
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

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-xmin")) { argv++; argc--; xscan_min = atoi(*argv); }
			else if (!strcmp(*argv, "-xmax")) { argv++; argc--; xscan_max = atoi(*argv); }
			else if (!strcmp(*argv, "-run")) { argv++; argc--; input_run = atoi(*argv); }
			else if (!strcmp(*argv, "-segment")) { argv++; argc--; input_segment = atoi(*argv); }
			else if (!strcmp(*argv, "-scan")) { argv++; argc--; input_scan = atoi(*argv); }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else if (!output_image_directory) output_image_directory = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
			argv++; argc--;
		}
	}

	// Check scene name
	if (!input_scene_name) {
		fprintf(stderr, "Usage: sadetect input_scene [output_directory] [-v] [-xmin] <xmin> [-xmax] <xmax> [-run] <run> [-segment] <segment> [-scan] <scan>\n");
		return FALSE;
	}

	// Check output image directory
	if (!output_image_directory) {
		output_image_directory = "gsv_data/laser_images";
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
