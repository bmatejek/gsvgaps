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
static char unaligned_scene_name[4096];
static char aligned_scene_name[4096];
static const char *input_run_name = NULL;
static const char *output_image_directory = NULL;
static int print_verbose = 0;
static char *object = NULL;

// Global variables
static GSVScene *unaligned_scene = NULL;
static GSVScene *aligned_scene = NULL;
static GSVRun *unaligned = NULL;
static GSVRun *aligned = NULL;
static vector<R3Point> unaligned_locations = vector<R3Point>();
static vector<R3Point> aligned_locations = vector<R3Point>();

// DA parameters
static RNLength DA_image_spacing = 0.1;

////////////////////////////////////////////////////////////////////////
// Grid names
////////////////////////////////////////////////////////////////////////

const char *object_name = "Truth";

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

static int DetectTruth(GSVScan *unaligned_scan, GSVScan *aligned_scan, const char *output_image_directory, const char *image_parameterization) {
	// Parameters
	RNScalar zoom_factor = 1;
	
	// Get convenient variables
	if (!unaligned_scan) return 1;
	if (unaligned_scan->NScanlines() == 0) return 1;
	int scan_index = unaligned_scan->SegmentIndex();

	// read in the points
	aligned_scan->ReadPoints();

	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = unaligned_scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;
	GSVScene *scene = run->Scene();
	const char *cache_directory = scene->CacheDataDirectoryName();
	char image_name[4096];

	
	// Print message
	if (print_verbose) {
		printf("    Creating original object truth images ...\n");
		fflush(stdout);
	}

	// read in unaligned position grids
	R2Grid unaligned_position_x_grid, unaligned_position_y_grid, unaligned_position_z_grid;
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_PositionX.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!unaligned_position_x_grid.Read(image_name)) return 0;
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_PositionY.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!unaligned_position_y_grid.Read(image_name)) return 0;
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_PositionZ.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!unaligned_position_z_grid.Read(image_name)) return 0;

	// read in the unaligned scanline position grids
	R2Grid unaligned_scanline_grid, unaligned_pointindex_grid;
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_Scanline.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!unaligned_scanline_grid.Read(image_name)) return 0;
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_PointIndex.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!unaligned_pointindex_grid.Read(image_name)) return 0;

	// Compute world to grid transformation
	int xres = (int)(zoom_factor * unaligned_scanline_grid.XResolution());
	int yres = (int)(zoom_factor * unaligned_scanline_grid.YResolution());

	// Create truth images for scan
	R2Affine unaligned_world_to_grid = R2identity_affine;
	unaligned_world_to_grid.Scale(zoom_factor);
	unaligned_world_to_grid.Transform(unaligned_scanline_grid.WorldToGridTransformation());
	R2Grid unaligned_truth_grid(xres, yres, unaligned_world_to_grid);
	unaligned_truth_grid.Clear(R2_GRID_UNKNOWN_VALUE);

	vector<R3Point> local_aligned_locations = vector<R3Point>();
	for (unsigned int t = 0; t < unaligned_locations.size(); t++) {
		R3Point pt = unaligned_locations[t];
		double d = DBL_MAX;
		int pxo = 0, pyo = 0;
		// go through every point to find the closest
		for (int px = 0; px < xres; px++) {
			for (int py = 0; py < yres; py++) {
				// get (x, y, z) world coordinates for this pixel
				RNScalar x = unaligned_position_x_grid.GridValue(px, py);
				if (x == R2_GRID_UNKNOWN_VALUE) continue;
				RNScalar y = unaligned_position_y_grid.GridValue(px, py);
				if (y == R2_GRID_UNKNOWN_VALUE) continue;
				RNScalar z = unaligned_position_z_grid.GridValue(px, py);
				if (z == R2_GRID_UNKNOWN_VALUE) continue;

				// see if this location is the closest seen so far
				R3Point pp = R3Point(x, y, z);
				if (R3Distance(pp, pt) < d) {
					d = R3Distance(pp, pt);
					pxo = px;
					pyo = py;
				}
			}
		}

		// make sure the distance is less than 1 meter away
		double epsilon = 1;
		if (d > epsilon) continue;

		// remove this location from the list of considerable locations
		unaligned_locations.erase(unaligned_locations.begin() + t); t--;

		// update this truth_grid
		unaligned_truth_grid.SetGridValue(pxo, pyo, 1);

		RNScalar scanline_index_value = unaligned_scanline_grid.GridValue(pxo, pyo);
		if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) { fprintf(stderr, "Failed to find corresponding scanline_grid value\n"); return 0; }
		RNScalar point_index_value = unaligned_pointindex_grid.GridValue(pxo, pyo);
		if (point_index_value == R2_GRID_UNKNOWN_VALUE) { fprintf(stderr, "Failed to find corresponding pointindex_grid value\n"); return 0; }

		int scanline_index = (int)(scanline_index_value + 0.5);
		int point_index = (int)(point_index_value + 0.5);

		GSVScanline *scanline = aligned_scan->Scanline(scanline_index);
		const R3Point& position = scanline->PointPosition(point_index);

		// local aligned ensures that no unaligned points correspond to multiple aligned points
		local_aligned_locations.push_back(R3Point(position));
		aligned_locations.push_back(R3Point(position));
	}

	// release the points
	aligned_scan->ReleasePoints();

	// Write scan grids 
	sprintf(image_name, "original/%s/laser_images/%s/%02d_%02d_%s_%s_Truth.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization, object);
	unaligned_truth_grid.Write(image_name);


	// Write aligned truth grid
	R2Grid aligned_position_x_grid, aligned_position_y_grid, aligned_position_z_grid;
	sprintf(image_name, "aligned/%s/laser_images/%s/%02d_%02d_%s_PositionX.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!aligned_position_x_grid.Read(image_name)) return 0;
	sprintf(image_name, "aligned/%s/laser_images/%s/%02d_%02d_%s_PositionY.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!aligned_position_y_grid.Read(image_name)) return 0;
	sprintf(image_name, "aligned/%s/laser_images/%s/%02d_%02d_%s_PositionZ.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!aligned_position_z_grid.Read(image_name)) return 0;

	// read in the aligned scanline position grids
	R2Grid aligned_scanline_grid, aligned_pointindex_grid;
	sprintf(image_name, "aligned/%s/laser_images/%s/%02d_%02d_%s_Scanline.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!aligned_scanline_grid.Read(image_name)) return 0;

	// Compute world to grid transformation
	xres = (int)(zoom_factor * aligned_scanline_grid.XResolution());
	yres = (int)(zoom_factor * aligned_scanline_grid.YResolution());

	// Create truth images for scan
	R2Affine aligned_world_to_grid = R2identity_affine;
	aligned_world_to_grid.Scale(zoom_factor);
	aligned_world_to_grid.Transform(aligned_scanline_grid.WorldToGridTransformation());
	R2Grid aligned_truth_grid(xres, yres, aligned_world_to_grid);
	aligned_truth_grid.Clear(R2_GRID_UNKNOWN_VALUE);

	for (unsigned int t = 0; t < local_aligned_locations.size(); t++) {
		R3Point pt = local_aligned_locations[t];
		double d = DBL_MAX;
		int pxo = 0, pyo = 0;
		// go through every point to find the closest
		for (int px = 0; px < xres; px++) {
			for (int py = 0; py < yres; py++) {
				// get (x, y, z) world coordinates for this pixel
				RNScalar x = aligned_position_x_grid.GridValue(px, py);
				if (x == R2_GRID_UNKNOWN_VALUE) continue;
				RNScalar y = aligned_position_y_grid.GridValue(px, py);
				if (y == R2_GRID_UNKNOWN_VALUE) continue;
				RNScalar z = aligned_position_z_grid.GridValue(px, py);
				if (z == R2_GRID_UNKNOWN_VALUE) continue;

				// see if this location is the closest seen so far
				R3Point pp = R3Point(x, y, z);
				if (R3Distance(pp, pt) < d) {
					d = R3Distance(pp, pt);
					pxo = px;
					pyo = py;
				}
			}
		}

		// update this truth_grid
		aligned_truth_grid.SetGridValue(pxo, pyo, 1);
	}

	// Write scan grids 
	sprintf(image_name, "aligned/%s/laser_images/%s/%02d_%02d_%s_%s_Truth.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization, object);
	aligned_truth_grid.Write(image_name);

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// SA images (Scanline vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteSAImages(GSVScan *unaligned_scan, GSVScan *aligned_scan, const char *output_image_directory) {
	// Get convenient variables
	if (!unaligned_scan) return 1;
	if (unaligned_scan->NScanlines() == 0) return 1;
	int scan_index = unaligned_scan->SegmentIndex();

	// don't care about the middle scan
	if (scan_index == 1) return 1;

	GSVSegment *segment = unaligned_scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Print message
	if (print_verbose) {
		printf("  Creating SA truth images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
		fflush(stdout);
	}

	// Capture SA image
	if (!DetectTruth(unaligned_scan, aligned_scan, output_image_directory, "SA")) return 0;

	// Return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// DA images (Distance vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteDAImage(GSVScan *scan, const char *output_image_directory, const R2Grid& sa_travel_distance_grid, const R2Grid& sa_viewpoint_distance_grid, const char *grid_name) {
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

	// get output file name
	char da_name[4096];
	sprintf(da_name, "%s/%s/%02d_%02d_DA_%s_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, object, grid_name);

	// Read SA image
	R2Grid sa_grid;
	char sa_name[4096];
	sprintf(sa_name, "%s/%s/%02d_%02d_SA_%s_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, object, grid_name);
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
	output_image_directory = "original/gsv_data/laser_images";
	sprintf(sa_travel_distance_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	sprintf(sa_viewpoint_distance_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);

	sa_travel_distance_grid.Read(sa_travel_distance_name);
	sa_viewpoint_distance_grid.Read(sa_viewpoint_distance_name);
	if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, object_name)) return 0;
	
	output_image_directory = "aligned/gsv_data/laser_images";
	sprintf(sa_travel_distance_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	sprintf(sa_viewpoint_distance_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	sa_travel_distance_grid.Read(sa_travel_distance_name);
	sa_viewpoint_distance_grid.Read(sa_viewpoint_distance_name);
	if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, object_name)) return 0;

	// Return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Top-level image computation function
////////////////////////////////////////////////////////////////////////

static int WriteImages(void) {
	// Start statistics
	RNTime start_time;
	start_time.Read();
	if (print_verbose) {
		printf("Creating images ...\n");
		fflush(stdout);
	}

	// Write scan images
	for (int ir = 0; ir < unaligned_scene->NRuns(); ir++) {
		GSVRun *unaligned_run = unaligned_scene->Run(ir);
		GSVRun *aligned_run = aligned_scene->Run(ir);
		for (int is = 0; is < unaligned_run->NSegments(); is++) {
			GSVSegment *unaligned_segment = unaligned_run->Segment(is);
			GSVSegment *aligned_segment = aligned_run->Segment(is);
			for (int ia = 0; ia < unaligned_segment->NScans(); ia++) {
				GSVScan *unaligned_scan = unaligned_segment->Scan(ia);
				GSVScan *aligned_scan = aligned_segment->Scan(ia);
				if (!WriteSAImages(unaligned_scan, aligned_scan, output_image_directory)) return 0;
				if (!WriteDAImages(unaligned_scan, output_image_directory)) return 0;
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
// Truth reading and writing functions
////////////////////////////////////////////////////////////////////////

static int ReadTruth(void) {
	char truthFilename[4096];
	sprintf(truthFilename, "original/truth/%s/%s/truth.txt", object, unaligned->Name());

	ifstream truth;
	truth.open(truthFilename);

	if (!truth.is_open()) { fprintf(stderr, "Failed to read unaligned truth: %s\n", truthFilename); return 0; }

	// read in all of the truth locations
	string line;
	while (getline(truth, line)) {
		vector<string> coordinates = splitString(line);
		if (coordinates.size() != 3) { fprintf(stderr, "Invalid coordinate in %s\n", truthFilename); return 0; }
		
		double x = atof(coordinates[0].c_str());
		double y = atof(coordinates[1].c_str());
		double z = atof(coordinates[2].c_str());

		unaligned_locations.push_back(R3Point(x, y, z));
	}

	truth.close();
	return 1;
}

static int SaveTruth(void) {
	char truthFilename[4096];
	sprintf(truthFilename, "aligned/truth/%s/%s/truth.txt", object, aligned->Name());

	ofstream truth;
	truth.open(truthFilename, ofstream::trunc);
	
	if (!truth.is_open()) { fprintf(stderr, "Failed to write aligned truth: %s\n", truthFilename);  return 0; }
	
	for (unsigned int it = 0; it < aligned_locations.size(); it++) {
		R3Point pt = aligned_locations[it];

		truth << pt.X();
		truth << ',';
		truth << pt.Y();
		truth << ',';
		truth << pt.Z();
		truth << '\n';
	}

	truth.close();
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage() {
	fprintf(stderr, "Usage: converttruth [run_name] -object <object>\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-object")) { argv++; argc--; object = *argv; }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_run_name) { input_run_name = *argv; }
			else if (!output_image_directory) output_image_directory = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// Check all input
	if (!input_run_name) return Usage();
	if (!object) return Usage();

	sprintf(unaligned_scene_name, "original/gsv_data/laser_points/%s/%s.gsv", input_run_name, input_run_name);
	sprintf(aligned_scene_name, "aligned/gsv_data/laser_points/%s/%s.gsv", input_run_name, input_run_name);



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

	// Read in scenes
	unaligned_scene = ReadScene(unaligned_scene_name);
	aligned_scene = ReadScene(aligned_scene_name);
	if (!unaligned_scene || !aligned_scene) exit(-1);

	// make sure scene are identical
	if (aligned_scene->NRuns() != unaligned_scene->NRuns() || unaligned_scene->NRuns() != 1) {
		fprintf(stderr, "Scenes must have one run\n");
		fprintf(stderr, "  Currently unaligned: %d\n", unaligned_scene->NRuns());
		fprintf(stderr, "  Currently aligned: %d\n", aligned_scene->NRuns());
		exit(-1);
	}

	// make sure correspond to the same scene name
	unaligned = unaligned_scene->Run(0);
	aligned = aligned_scene->Run(0);

	// read in unaligned truth
	if (!ReadTruth()) exit(-1);

	// create truth grids
	unsigned int unaligned_locations_size = unaligned_locations.size();
	if (!WriteImages()) exit(-1);

	// confirm that no points were lose
	if (unaligned_locations_size != aligned_locations.size()) {
		fprintf(stderr, "%s %s\n", unaligned->Name(), object);
		fprintf(stderr, "Points lost in transition\n");
		fprintf(stderr, "  Unaligned true locations: %u\n", unaligned_locations_size);
		fprintf(stderr, "  Aligned true locations: %lu\n", aligned_locations.size());
	}

	// save the aligned locations
	if (!SaveTruth()) exit(-1);

	// Return success 
	return 0;
}