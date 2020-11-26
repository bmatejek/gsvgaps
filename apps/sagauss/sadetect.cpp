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
static char *objects_file_name = "rcnn/objects.txt";
static int print_verbose = 0;

// Depth and color parameters
static double min_frustum_depth = 1;
static double max_frustum_depth = 1000;

// Global variables
static int GSV_IMAGE_WIDTH = 1936;
static int GSV_IMAGE_HEIGHT = 2592;
static double low_threshold = 0;
static double high_threshold = 1.0;
map<string, int> object_map;			// hash function to map object strings to ids

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
	int xres = (int)(detection_zoom_factor * scanline_grid.XResolution());
	int yres = (int)(detection_zoom_factor * scanline_grid.YResolution());

	// Create detection images for scan
	R2Affine world_to_grid = R2identity_affine;
	world_to_grid.Scale(detection_zoom_factor);
	world_to_grid.Transform(scanline_grid.WorldToGridTransformation());

	// go through every image
	for (int ir = 0; ir < scene->NRuns(); ++ir) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0, pn = 0; is < run->NSegments(); ++is) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++pn) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ++ii) {
					printf("%d_%d\n", ip, ii);
					GSVImage *image = panorama->Image(ii);

					// open up RCNN file
					char RCNNFilename[4096];
					sprintf(RCNNFilename, "rcnn/%s/%02d_%06d_%02d_UndistortedImage.txt", run->Name(), is, pn, ii);
					ifstream RCNNFile;
					RCNNFile.open(RCNNFilename);
					if (!RCNNFile.is_open()) { fprintf(stderr, "Failed to read %s\n", RCNNFilename); return 0; }

					// set id to index mapping
					map<int, int> idMap = map<int, int>();	
					string line;

					while (getline(RCNNFile, line)) {
						vector<string> box = splitString(line);

						string object = box[0];
						int id = object_map[object];
						idMap.insert(pair<int, int>())
					}

					RCNNFile.clear();
					
					// read into a grid
					R2Grid **values = new R2Grid *[idMap.size()];

					// create new R2Grids
					for (int i = 0; i < idMap.size(); i++) {
						values[i] = new R2Grid(GSV_IMAGE_WIDTH, GSV_IMAGE_HEIGHT);
						for (int x = 0; x < values[i]->XResolution(); x++) {
							for (int y = 0; y < values[i]->YResolution(); y++) {
								values[i]->SetGridValue(x, y, R2_GRID_UNKNOWN_VALUE);
							}
						}
					}

					while (getline(RCNNFile, line)) {
						vector<string> box = splitString(line);

						string object = box[0];
						int id = object_map[object];
						double score = atof(box[1].c_str());
						double xmin = atof(box[2].c_str());
						double ymax = GSV_IMAGE_HEIGHT - atof(box[3].c_str());
						double xmax = atof(box[4].c_str());
						double ymin = GSV_IMAGE_HEIGHT - atof(box[5].c_str());

						int index = idMap[id];

						for (int x = xmin; x < xmax; x++) {
							for (int y = ymin; y < ymax; y++) {
								if (score > values[index]->GridValue(x, y)) {
									values[index]->SetGridValue(x, y, score);
								}
							}
						}
					}

					R2Grid detection_grid(xres, yres, world_to_grid);
					detection_grid.Clear(R2_GRID_UNKNOWN_VALUE);

					// project to each point
					for (int ix = 0; ix < xres; ix++) {
						for (int iy = 0; iy < yres; iy++) {
							// find this pixel
							int px = (int)(ix / detection_zoom_factor + 0.5); if (px >= position_x_grid.XResolution()) continue;
							int py = (int)(iy / detection_zoom_factor + 0.5); if (py >= position_x_grid.YResolution()) continue;

							// get (x, y, z) world coordinates of this pixel
							RNScalar x = position_x_grid.GridValue(px, py); if (x == R2_GRID_UNKNOWN_VALUE) continue;
							RNScalar y = position_y_grid.GridValue(px, py);	if (y == R2_GRID_UNKNOWN_VALUE) continue;
							RNScalar z = position_z_grid.GridValue(px, py);	if (z == R2_GRID_UNKNOWN_VALUE) continue;

							// make sure there is a scanline here
							RNScalar scanline_value = scanline_grid.GridValue(px, py); 
							if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
							int scanline_index = (int)(scanline_value + 0.5);

							// get scanline
							GSVScanline *scanline = scan->Scanline(scanline_index);
							R3Point position(x, y, z);

							// Get image info
							const GSVCamera *camera = image->Camera();
							const GSVPose& image_pose = image->Pose();
							const R3Point& image_viewpoint = image_pose.Viewpoint();

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


						}
					}

					// free up memory
					delete values;

					// close the file, no longer needed
					RCNNFile.close();
				}
			}
		}
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
		fprintf(stderr, "Usage: sagauss input_scene [output_directory] [-v]\n");
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

	// Read in all object file
	if (!ReadObjectFile()) { exit(-1); }

	// Write images
	if (!WriteImages(scene, output_image_directory)) exit(-1);

	// Return success 
	return 0;
}