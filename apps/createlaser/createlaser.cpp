////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments
static const char *input_scene_name = NULL;
static int print_verbose = 0;
static int overwrite_all = 0;

// global variables
char *output_image_directory = (char *) "gsv_data/laser_images";

// DA parameters
static RNLength SA_viewpoint_spacing = 0.01;
static RNLength DA_image_spacing = 0.1;

////////////////////////////////////////////////////////////////////////
// Grid names
////////////////////////////////////////////////////////////////////////

// Grid names
const char *base_image_names[] = {
	"Scanline", "ViewpointDistance", "PositionX", "PositionY", "PositionZ", "PointIndex", "Timestamp"
};

int num_base_image_names = sizeof(base_image_names) / sizeof(const char *);

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

////////////////////////////////////////////////////////////////////////
// SA images (Scanline vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteSABaseImages(GSVScan *scan, const char *output_image_directory) {

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

	// Skip if already done
	char image_name[4096];
	sprintf(image_name, "%s/%s/%02d_%02d_SA_Timestamp.grd", output_image_directory, run->Name(), segment_index, scan_index);
	if (RNFileExists(image_name) && !overwrite_all) return 1;

	// Print message
	if (print_verbose) {
		printf("    Creating base images ...\n");
		fflush(stdout);
	}

	// Count scanlines
	int nscanlines = 0;
	R3Point prev_viewpoint = scan->Scanline(0)->Pose().Viewpoint();
	for (int ie = 0; ie < scan->NScanlines(); ie++) {
		GSVScanline *scanline = scan->Scanline(ie);
		const GSVPose& pose = scanline->Pose();
		const R3Point& viewpoint = pose.Viewpoint();
		if ((ie > 0) && (R3Distance(viewpoint, prev_viewpoint) < SA_viewpoint_spacing)) continue;
		prev_viewpoint = viewpoint;
		nscanlines++;
	}

	// Compute grid resolution
	int xres = nscanlines;
	int yres = 180;

	// Initialize grids
	R2Box grid_box(0, 0, scan->NScanlines(), 180);
	R2Grid timestamp_grid(xres, yres, grid_box);       timestamp_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid distance_grid(xres, yres, grid_box);        distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid travel_distance_grid(xres, yres, grid_box); travel_distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid scanline_grid(xres, yres, grid_box);        scanline_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid pointindex_grid(xres, yres, grid_box);      pointindex_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionX_grid(xres, yres, grid_box);       positionX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionY_grid(xres, yres, grid_box);       positionY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionZ_grid(xres, yres, grid_box);       positionZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	
	// Read scan points
	if (!scan->ReadPoints()) {
		fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
		return 0;
	}

	// Fill grids
	nscanlines = 0;
	RNScalar travel_distance = 0;
	prev_viewpoint = scan->Scanline(0)->Pose().Viewpoint();
	for (int ie = 0; ie < scan->NScanlines(); ie++) {
		GSVScanline *scanline = scan->Scanline(ie);
		const GSVPose& pose = scanline->Pose();
		const R3Point& viewpoint = pose.Viewpoint();
		if ((ie > 0) && (R3Distance(viewpoint, prev_viewpoint) < SA_viewpoint_spacing)) continue;
		travel_distance += R3Distance(viewpoint, prev_viewpoint);
		prev_viewpoint = viewpoint;
		int ix = nscanlines++;
		if (ix >= xres) continue;
		RNScalar timestamp = scanline->Timestamp();
		const R3Vector& towards = pose.Towards();
		const R3Vector& up = pose.Up();
		for (int j = 0; j < scanline->NPoints(); j++) {
			R3Point position = scanline->PointPosition(j);
			R3Vector v = position - viewpoint;
			RNLength distance = v.Length();
			if (RNIsZero(distance)) continue;
			v /= distance;
			RNScalar dot = v.Dot(towards);
			RNAngle angle = (dot < 1) ? ((dot > -1) ? acos(dot) : RN_PI) : 0;
			if (v.Dot(up) > 0) angle = RN_PI_OVER_TWO + angle;
			else angle = RN_PI_OVER_TWO - angle;
			int iy = (int)(yres * angle / RN_PI + 0.5);
			if ((iy < 0) || (iy >= yres)) continue;

			// Remember values for closest scan point
			RNScalar old_distance_value = distance_grid.GridValue(ix, iy);
			if ((old_distance_value == R2_GRID_UNKNOWN_VALUE) || (distance < old_distance_value)) {
				timestamp_grid.SetGridValue(ix, iy, timestamp);
				distance_grid.SetGridValue(ix, iy, distance);
				travel_distance_grid.SetGridValue(ix, iy, travel_distance);
				scanline_grid.SetGridValue(ix, iy, ie);
				pointindex_grid.SetGridValue(ix, iy, j);
				positionX_grid.SetGridValue(ix, iy, position.X());
				positionY_grid.SetGridValue(ix, iy, position.Y());
				positionZ_grid.SetGridValue(ix, iy, position.Z());
			}
		}
	}

	// Release scan points
	if (!scan->ReleasePoints()) {
		fprintf(stderr, "Unable to release points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
		return 0;
	}

	// Create hole grid
	R2Grid hole_grid(distance_grid);
	hole_grid.Clear(0);
	hole_grid.Substitute(R2_GRID_UNKNOWN_VALUE, 1);

	// Write grids
	sprintf(image_name, "%s/%s/%02d_%02d_SA_Timestamp.grd", output_image_directory, run->Name(), segment_index, scan_index);
	timestamp_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	distance_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
	travel_distance_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_Scanline.grd", output_image_directory, run->Name(), segment_index, scan_index);
	scanline_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_PointIndex.grd", output_image_directory, run->Name(), segment_index, scan_index);
	pointindex_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index);
	positionX_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index);
	positionY_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
	positionZ_grid.Write(image_name);
	
	// Return success
	return 1;
}

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
	if (!WriteSABaseImages(scan, output_image_directory)) return 0;

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

	// Skip if already done
	char da_name[4096];
	sprintf(da_name, "%s/%s/%02d_%02d_DA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
	if (RNFileExists(da_name) && !overwrite_all) return 1;

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
	for (int i = 0; i < num_base_image_names; i++) {
		if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, base_image_names[i])) return 0;
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
	fprintf(stderr, "Usage: createlaser input_scene -v [-overwrite_all]\n");
	return FALSE;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-overwrite_all")) { overwrite_all = 1;  }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// Check all input
	if (!input_scene_name) return Usage();

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

	// Write images
	if (!WriteImages(scene, output_image_directory)) exit(-1);

	// Return success 
	return 0;
}