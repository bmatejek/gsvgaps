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
static int overwrite_detection = 0;

// XY, DA, DH parameters
static RNLength SA_viewpoint_spacing = 0.01;

// Depth and color parameters
static double min_frustum_depth = 1;
static double max_frustum_depth = 1000;
static double max_viewpoint_distance = 100;
static double max_timestamp_difference = 100;

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
// Parameterization independent grid functions
////////////////////////////////////////////////////////////////////////

static int WriteNormalImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization) {
	// Get convenient variables
	if (!scan) return 0;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Skip if already done
	char image_name[4096];
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (RNFileExists(image_name)) return 1;

	// Print message
	if (print_verbose) {
		printf("    Creating normal images ...\n");
		fflush(stdout);
	}

	// Read position grids
	R2Grid positionX_grid, positionY_grid, positionZ_grid;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionX_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionY_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionZ_grid.Read(image_name)) return 0;

	// Read viewpoint grids
	R2Grid viewpointX_grid, viewpointY_grid, viewpointZ_grid;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!viewpointX_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!viewpointY_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!viewpointZ_grid.Read(image_name)) return 0;

	// Create normal grids
	R2Grid normalX_grid(positionX_grid);   normalX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid normalY_grid(normalX_grid);
	R2Grid normalZ_grid(normalX_grid);

	// Fill normal grids
	int xres = normalX_grid.XResolution();
	int yres = normalX_grid.YResolution();
	for (int i = 0; i < xres; i++) {
		int i0 = (i > 0) ? i - 1 : 0;
		int i1 = (i < xres - 1) ? i + 1 : xres - 1;
		for (int j = 0; j < yres; j++) {
			int j0 = (j > 0) ? j - 1 : 0;
			int j1 = (j < yres - 1) ? j + 1 : yres - 1;

			// Check position values
			if (positionX_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
			if (positionY_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
			if (positionZ_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;

			// Get 3D vector between horizontal neighbors
			RNScalar xA0 = positionX_grid.GridValue(i0, j);
			if (xA0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar yA0 = positionY_grid.GridValue(i0, j);
			if (yA0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar zA0 = positionZ_grid.GridValue(i0, j);
			if (zA0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar xA1 = positionX_grid.GridValue(i1, j);
			if (xA1 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar yA1 = positionY_grid.GridValue(i1, j);
			if (yA1 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar zA1 = positionZ_grid.GridValue(i1, j);
			if (zA1 == R2_GRID_UNKNOWN_VALUE) continue;
			R3Vector vA(xA1 - xA0, yA1 - yA0, zA1 - zA0);
			vA.Normalize();

			// Get 3D vector between vertical neighbors
			RNScalar xB0 = positionX_grid.GridValue(i, j0);
			if (xB0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar yB0 = positionY_grid.GridValue(i, j0);
			if (yB0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar zB0 = positionZ_grid.GridValue(i, j0);
			if (zB0 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar xB1 = positionX_grid.GridValue(i, j1);
			if (xB1 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar yB1 = positionY_grid.GridValue(i, j1);
			if (yB1 == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar zB1 = positionZ_grid.GridValue(i, j1);
			if (zB1 == R2_GRID_UNKNOWN_VALUE) continue;
			R3Vector vB(xB1 - xB0, yB1 - yB0, zB1 - zB0);
			vB.Normalize();

			// Compute view vector
			RNScalar dx = positionX_grid.GridValue(i, j) - viewpointX_grid.GridValue(i, j);
			RNScalar dy = positionY_grid.GridValue(i, j) - viewpointY_grid.GridValue(i, j);
			RNScalar dz = positionZ_grid.GridValue(i, j) - viewpointZ_grid.GridValue(i, j);
			R3Vector view_vector(dx, dy, dz);
			view_vector.Normalize();

			// Compute normal
			R3Vector normal = vA % vB;
			normal.Normalize();

			// Update normal and NdotV grids
			normalX_grid.SetGridValue(i, j, normal.X());
			normalY_grid.SetGridValue(i, j, normal.Y());
			normalZ_grid.SetGridValue(i, j, normal.Z());
		}
	}

	// Write grids
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	normalX_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	normalY_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	normalZ_grid.Write(image_name);

	// Return success
	return 1;
}

static int WriteCurvatureImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization) {
	// Get convenient variables
	if (!scan) return 0;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Skip if already done
	char image_name[4096];
	sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (RNFileExists(image_name)) return 1;

	// Print message
	if (print_verbose) {
		printf("    Creating curvature images ...\n");
		fflush(stdout);
	}

	// Read position grids
	R2Grid positionX_grid, positionY_grid, positionZ_grid;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionX_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionY_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!positionZ_grid.Read(image_name)) return 0;
	positionX_grid.Blur(1);
	positionY_grid.Blur(1);
	positionZ_grid.Blur(1);

	// Read normal grids
	R2Grid normalX_grid, normalY_grid, normalZ_grid;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!normalX_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!normalY_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!normalZ_grid.Read(image_name)) return 0;
	normalX_grid.Blur(1);
	normalY_grid.Blur(1);
	normalZ_grid.Blur(1);

	// Create curvature grids
	R2Grid horizontal_curvature_grid(positionX_grid);
	horizontal_curvature_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid vertical_curvature_grid(horizontal_curvature_grid);

	// Compute whether left-right directions should be flipped 
	int dir = (scan->SegmentIndex() != 2) ? 1 : -1;

	// Fill horizontal curvature grid
	for (int i = 1; i < horizontal_curvature_grid.XResolution() - 1; i++) {
		for (int j = 0; j < horizontal_curvature_grid.YResolution(); j++) {
			// Get position values
			RNScalar positionX = positionX_grid.GridValue(i, j);
			if (positionX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionY = positionY_grid.GridValue(i, j);
			if (positionY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionZ = positionZ_grid.GridValue(i, j);
			if (positionZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Point position(positionX, positionY, positionZ);

			// Get normal values
			RNScalar normalX = normalX_grid.GridValue(i, j);
			if (normalX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar normalY = normalY_grid.GridValue(i, j);
			if (normalY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar normalZ = normalZ_grid.GridValue(i, j);
			if (normalZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Vector normal(normalX, normalY, normalZ);

			// Get neighbor positions
			RNScalar positionAX = positionX_grid.GridValue(i - dir, j);
			if (positionAX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionAY = positionY_grid.GridValue(i - dir, j);
			if (positionAY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionAZ = positionZ_grid.GridValue(i - dir, j);
			if (positionAZ == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBX = positionX_grid.GridValue(i + dir, j);
			if (positionBX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBY = positionY_grid.GridValue(i + dir, j);
			if (positionBY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBZ = positionZ_grid.GridValue(i + dir, j);
			if (positionBZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Point positionA(positionAX, positionAY, positionAZ);
			R3Point positionB(positionBX, positionBY, positionBZ);
			R3Vector vA = positionA - position;
			R3Vector vB = positionB - position;
			RNLength lenA = vA.Length();
			RNLength lenB = vB.Length();
			if (RNIsZero(lenA) || RNIsZero(lenB)) continue;
			vA /= lenA; vB /= lenB;
			if (lenA < 1) lenA = 1;
			if (lenB < 1) lenB = 1;
			RNScalar angle = R3InteriorAngle(vA, vB);
			R3Vector dir_vector = vB + vA;
			RNScalar dot = dir_vector.Dot(normal);
			RNScalar sign = (dot > 0) ? -1 : 1;
			RNScalar curvature = sign * (RN_PI - angle) / (lenA + lenB);
			horizontal_curvature_grid.SetGridValue(i, j, curvature);
		}
	}

	// Fill vertical curvature grid
	for (int i = 0; i < horizontal_curvature_grid.XResolution(); i++) {
		for (int j = 1; j < horizontal_curvature_grid.YResolution() - 1; j++) {
			// Get position values
			RNScalar positionX = positionX_grid.GridValue(i, j);
			if (positionX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionY = positionY_grid.GridValue(i, j);
			if (positionY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionZ = positionZ_grid.GridValue(i, j);
			if (positionZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Point position(positionX, positionY, positionZ);

			// Get normal values
			RNScalar normalX = normalX_grid.GridValue(i, j);
			if (normalX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar normalY = normalY_grid.GridValue(i, j);
			if (normalY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar normalZ = normalZ_grid.GridValue(i, j);
			if (normalZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Vector normal(normalX, normalY, normalZ);

			// Get neighbor positions
			RNScalar positionAX = positionX_grid.GridValue(i, j - 1);
			if (positionAX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionAY = positionY_grid.GridValue(i, j - 1);
			if (positionAY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionAZ = positionZ_grid.GridValue(i, j - 1);
			if (positionAZ == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBX = positionX_grid.GridValue(i, j + 1);
			if (positionBX == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBY = positionY_grid.GridValue(i, j + 1);
			if (positionBY == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar positionBZ = positionZ_grid.GridValue(i, j + 1);
			if (positionBZ == R2_GRID_UNKNOWN_VALUE) continue;
			R3Point positionA(positionAX, positionAY, positionAZ);
			R3Point positionB(positionBX, positionBY, positionBZ);
			R3Vector vA = positionA - position;
			R3Vector vB = positionB - position;
			RNLength lenA = vA.Length();
			RNLength lenB = vB.Length();
			if (RNIsZero(lenA) || RNIsZero(lenB)) continue;
			vA /= lenA; vB /= lenB;
			if (lenA < 1) lenA = 1;
			if (lenB < 1) lenB = 1;
			RNScalar angle = R3InteriorAngle(vA, vB);
			R3Vector dir_vector = vB + vA;
			RNScalar dot = dir_vector.Dot(normal);
			RNScalar sign = (dot > 0) ? -1 : 1;
			RNScalar curvature = sign * (RN_PI - angle) / (lenA + lenB);
			vertical_curvature_grid.SetGridValue(i, j, curvature);
		}
	}

	// Write grids
	sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	horizontal_curvature_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_VerticalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	vertical_curvature_grid.Write(image_name);

	// Return success
	return 1;
}

static int WriteBoundaryImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization, RNLength max_depth_discontinuity = 1.0, RNScalar curvature_threshold = 0.25) {
	// Get convenient variables
	if (!scan) return 0;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Skip if already done
	char image_name[4096];
	sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (RNFileExists(image_name)) return 1;

	// Print message
	if (print_verbose) {
		printf("    Creating boundary images ...\n");
		fflush(stdout);
	}

	// Compute whether left-right directions should be flipped 
	int dir = (scan->SegmentIndex() != 2) ? 1 : -1;

	// Read viewpoint depth grid
	R2Grid viewpoint_depth_grid;
	R2Grid horizontal_curvature_grid;
	R2Grid vertical_curvature_grid;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointDepth.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!viewpoint_depth_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!horizontal_curvature_grid.Read(image_name)) return 0;
	sprintf(image_name, "%s/%s/%02d_%02d_%s_VerticalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (!vertical_curvature_grid.Read(image_name)) return 0;

	// Create boundary grid
	R2Grid boundary_type_grid(viewpoint_depth_grid.XResolution(), viewpoint_depth_grid.YResolution());

	// Fill boundary grids
	for (int i = 1; i < boundary_type_grid.XResolution() - 1; i++) {
		for (int j = 1; j < boundary_type_grid.YResolution() - 1; j++) {
			RNScalar depth = viewpoint_depth_grid.GridValue(i, j);
			if (depth == R2_GRID_UNKNOWN_VALUE) continue;

			// Initialize boundary type 
			int boundary_type = 0;
			RNScalar neighbor_depth;

			// Check left
			neighbor_depth = viewpoint_depth_grid.GridValue(i - dir, j);
			if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_LEFT_UNKNOWN_BOUNDARY | GSV_LEFT_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_LEFT_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_LEFT_SHADOW_BOUNDARY;

			// Check right
			neighbor_depth = viewpoint_depth_grid.GridValue(i + dir, j);
			if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_RIGHT_UNKNOWN_BOUNDARY | GSV_RIGHT_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_RIGHT_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_RIGHT_SHADOW_BOUNDARY;

			// Check down
			neighbor_depth = viewpoint_depth_grid.GridValue(i, j - 1);
			if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_DOWN_UNKNOWN_BOUNDARY | GSV_DOWN_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_DOWN_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_DOWN_SHADOW_BOUNDARY;

			// Check up
			neighbor_depth = viewpoint_depth_grid.GridValue(i, j + 1);
			if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_UP_UNKNOWN_BOUNDARY | GSV_UP_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_UP_SILHOUETTE_BOUNDARY;
			else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_UP_SHADOW_BOUNDARY;

			// Check other types
			if (boundary_type == 0) {
				// Check horizontal curvature
				RNScalar hcurv = horizontal_curvature_grid.GridValue(i, j);
				if (hcurv != R2_GRID_UNKNOWN_VALUE) {
					if (hcurv > curvature_threshold) boundary_type += GSV_HORIZONTAL_RIDGE_BOUNDARY;
					else if (hcurv < -curvature_threshold) boundary_type += GSV_HORIZONTAL_VALLEY_BOUNDARY;
				}

				// Check vertical curvature
				RNScalar vcurv = vertical_curvature_grid.GridValue(i, j);
				if (vcurv != R2_GRID_UNKNOWN_VALUE) {
					if (vcurv > curvature_threshold) boundary_type += GSV_VERTICAL_RIDGE_BOUNDARY;
					else if (vcurv < -curvature_threshold) boundary_type += GSV_VERTICAL_VALLEY_BOUNDARY;
				}
			}

			// Set grid value
			boundary_type_grid.SetGridValue(i, j, boundary_type);
		}
	}

	// Write grids
	sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	boundary_type_grid.Write(image_name);

	// Return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Parameterization independent grid functions
////////////////////////////////////////////////////////////////////////

static int WriteColorImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization) {
	// Parameters
	RNScalar color_image_zoom_factor = 1;

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

	// Skip if already done
	sprintf(image_name, "%s/%s/%02d_%02d_%s_Color.bmp",
		output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	if (RNFileExists(image_name)) return 1;

	// Print message
	if (print_verbose) {
		printf("    Creating color images ...\n");
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
	int xres = (int)(color_image_zoom_factor * scanline_grid.XResolution());
	int yres = (int)(color_image_zoom_factor * scanline_grid.YResolution());

	// Get mesh (commented out just for efficiency)
	// GSVMesh *mesh = scan->Mesh();
	GSVMesh *mesh = NULL;

	// Create color images for scan
	R2Affine world_to_grid = R2identity_affine;
	world_to_grid.Scale(color_image_zoom_factor);
	world_to_grid.Transform(scanline_grid.WorldToGridTransformation());
	R2Grid scan_red_grid(xres, yres, world_to_grid);
	scan_red_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid scan_green_grid(scan_red_grid);
	R2Grid scan_blue_grid(scan_red_grid);
	R2Grid scan_timestamp_grid(scan_red_grid);

	// Consider tapestries relevant to scan
	int ntapestries[3] = { 3, 3, 3 };
	int tapestries[3][3] = { { 5, 6, 7 }, { 0, 1, 7 }, { 1, 2, 3 } };
	for (int k = 0; k < ntapestries[scan_index]; k++) {
		int it = tapestries[scan_index][k];
		GSVTapestry *tapestry = segment->Tapestry(it);
		GSVCamera *camera = tapestry->Camera();
		if (!camera) continue;

		// Print message
		if (print_verbose) {
			printf("      Creating color images for tapestry %02d\n", it);
			fflush(stdout);
		}

		// Consider every image in tapestry
		for (int ii = 0; ii < tapestry->NImages(); ii++) {
			GSVImage *image = tapestry->Image(ii);

			// Get image info
			const GSVPose& image_pose = image->Pose();
			RNScalar image_timestamp = image->Timestamp();
			const R3Point& image_viewpoint = image_pose.Viewpoint();
			R3Frustum image_frustum(image_viewpoint, image_pose.Towards(), image_pose.Up(),
				camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);

			// Get distorted image
			R2Image *distorted_image = image->DistortedImage();
			if (!distorted_image) continue;

			// Consider every grid pixel
			for (int ix = 0; ix < xres; ix++) {
				for (int iy = 0; iy < yres; iy++) {
					int px = (int)(ix / color_image_zoom_factor + 0.5);
					if (px >= position_x_grid.XResolution()) continue;
					int py = (int)(iy / color_image_zoom_factor + 0.5);
					if (py >= position_x_grid.YResolution()) continue;
					RNScalar x = position_x_grid.GridValue(px, py);
					if (x == R2_GRID_UNKNOWN_VALUE) continue;
					RNScalar y = position_y_grid.GridValue(px, py);
					if (y == R2_GRID_UNKNOWN_VALUE) continue;
					RNScalar z = position_z_grid.GridValue(px, py);
					if (z == R2_GRID_UNKNOWN_VALUE) continue;
					RNScalar scanline_value = scanline_grid.GridValue(px, py);
					if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
					int scanline_index = (int)(scanline_value + 0.5);
					GSVScanline *scanline = scan->Scanline(scanline_index);
					R3Point world_position(x, y, z);

					// Check if scanline timestamp is within range
					RNScalar scanline_timestamp = scanline->Timestamp();
					if (max_timestamp_difference > 0) {
						RNScalar timestamp_difference = fabs(scanline_timestamp - image->Timestamp());
						if (timestamp_difference > max_timestamp_difference) continue;
					}

					// Check if scanline viewpoint is within range
					const GSVPose& scanline_pose = scanline->Pose();
					RNLength viewpoint_distance = R3Distance(scanline_pose.Viewpoint(), image_viewpoint);
					if ((max_viewpoint_distance > 0) && (viewpoint_distance > max_viewpoint_distance)) continue;

					// Check if point is within image view frustum
					if (!image_frustum.Intersects(world_position)) continue;

					// Compute whether point is occluded
					if (mesh) {
						R3MeshIntersection intersection;
						R3Ray ray(image_viewpoint, world_position);
						const RNLength visibility_tolerance = 1.0;
						// if (mesh->Intersection(ray, &intersection, 0, viewpoint_distance)) {
						if (mesh->Intersection(ray, &intersection)) {
							if (intersection.t < viewpoint_distance - visibility_tolerance) continue;
						}
					}

					// Retrieve color of point from distorted image
					R2Point image_position = image->DistortedPosition(world_position);
					if (image_position.X() < 0) continue;
					if (image_position.Y() < 0) continue;
					int image_ix1 = (int)image_position.X();
					int image_iy1 = (int)image_position.Y();
					if ((image_ix1 < 0) || (image_ix1 >= distorted_image->Width())) continue;
					if ((image_iy1 < 0) || (image_iy1 >= distorted_image->Height())) continue;
					int image_ix2 = image_ix1 + 1;
					int image_iy2 = image_iy1 + 1;
					if ((image_ix2 < 0) || (image_ix2 >= distorted_image->Width())) continue;
					if ((image_iy2 < 0) || (image_iy2 >= distorted_image->Height())) continue;
					RNRgb color11 = distorted_image->PixelRGB(image_ix1, image_iy1);
					RNRgb color12 = distorted_image->PixelRGB(image_ix1, image_iy2);
					RNRgb color21 = distorted_image->PixelRGB(image_ix2, image_iy1);
					RNRgb color22 = distorted_image->PixelRGB(image_ix2, image_iy2);
					RNScalar tx = image_position.X() - image_ix1;
					RNScalar ty = image_position.Y() - image_iy1;
					RNRgb colorA = (1 - tx)*color11 + tx*color21;
					RNRgb colorB = (1 - tx)*color12 + tx*color22;
					RNRgb color = (1 - ty)*colorA + ty*colorB;

					// Update scan grid
					RNScalar scan_timestamp = scan_timestamp_grid.GridValue(ix, iy);
					if ((scan_timestamp == R2_GRID_UNKNOWN_VALUE) ||
						(fabs(scanline_timestamp - image_timestamp) < fabs(scanline_timestamp - scan_timestamp))) {
						scan_red_grid.SetGridValue(ix, iy, color.R());
						scan_green_grid.SetGridValue(ix, iy, color.G());
						scan_blue_grid.SetGridValue(ix, iy, color.B());
						scan_timestamp_grid.SetGridValue(ix, iy, image_timestamp);
					}
				}
			}

			// Delete distorted image
			delete distorted_image;
		}
	}

	// Create scan color image
	R2Image scan_rgb_image(xres, yres);
	for (int ix = 0; ix < xres; ix++) {
		for (int iy = 0; iy < yres; iy++) {
			RNScalar r = scan_red_grid.GridValue(ix, iy);
			if (r == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar g = scan_green_grid.GridValue(ix, iy);
			if (g == R2_GRID_UNKNOWN_VALUE) continue;
			RNScalar b = scan_blue_grid.GridValue(ix, iy);
			if (b == R2_GRID_UNKNOWN_VALUE) continue;
			RNRgb color(r, g, b);
			scan_rgb_image.SetPixelRGB(ix, iy, color);
		}
	}

	// Write scan grids 
	sprintf(image_name, "%s/%s/%02d_%02d_%s_Color.bmp",
		output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	scan_rgb_image.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorRed.grd",
		output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	scan_red_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorGreen.grd",
		output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	scan_green_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorBlue.grd",
		output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
	scan_blue_grid.Write(image_name);

	// Delete mesh
	if (mesh) delete mesh;

	// Return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// SA images (Scanline vs. Angle)
////////////////////////////////////////////////////////////////////////

static int WriteSABaseImages(GSVScan *scan, const char *output_image_directory) {
	// Get convenient variables
	if (!scan) return 1;
	if (scan->NScanlines() == 0) return 1;
	int scan_index = scan->SegmentIndex();
	GSVSegment *segment = scan->Segment();
	if (!segment) return 0;
	int segment_index = segment->RunIndex();
	GSVRun *run = segment->Run();
	if (!run) return 0;

	// Skip if already done
	char image_name[4096];
	sprintf(image_name, "%s/%s/%02d_%02d_SA_Timestamp.grd", output_image_directory, run->Name(), segment_index, scan_index);
	if (RNFileExists(image_name)) return 1;

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
	R2Grid depth_grid(xres, yres, grid_box);           depth_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid distance_grid(xres, yres, grid_box);        distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid timestamp_grid(xres, yres, grid_box);       timestamp_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid scanline_grid(xres, yres, grid_box);        scanline_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid pointindex_grid(xres, yres, grid_box);      pointindex_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionX_grid(xres, yres, grid_box);       positionX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionY_grid(xres, yres, grid_box);       positionY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid positionZ_grid(xres, yres, grid_box);       positionZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid viewpointX_grid(xres, yres, grid_box);      viewpointX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid viewpointY_grid(xres, yres, grid_box);      viewpointY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
	R2Grid viewpointZ_grid(xres, yres, grid_box);      viewpointZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);

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
			RNScalar depth = v.Dot(towards);
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
				depth_grid.SetGridValue(ix, iy, depth);
				distance_grid.SetGridValue(ix, iy, distance);
				scanline_grid.SetGridValue(ix, iy, ie);
				pointindex_grid.SetGridValue(ix, iy, j);
				positionX_grid.SetGridValue(ix, iy, position.X());
				positionY_grid.SetGridValue(ix, iy, position.Y());
				positionZ_grid.SetGridValue(ix, iy, position.Z());
				viewpointX_grid.SetGridValue(ix, iy, viewpoint.X());
				viewpointY_grid.SetGridValue(ix, iy, viewpoint.Y());
				viewpointZ_grid.SetGridValue(ix, iy, viewpoint.Z());
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
	sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointDepth.grd", output_image_directory, run->Name(), segment_index, scan_index);
	depth_grid.Write(image_name);
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
	sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointX.grd", output_image_directory, run->Name(), segment_index, scan_index);
	viewpointX_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointY.grd", output_image_directory, run->Name(), segment_index, scan_index);
	viewpointY_grid.Write(image_name);
	sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
	viewpointZ_grid.Write(image_name);

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
	if (!WriteSABaseImages(scan, output_image_directory)) return 0;
	if (!WriteNormalImages(scan, output_image_directory, "SA")) return 0;
	if (!WriteCurvatureImages(scan, output_image_directory, "SA")) return 0;
	if (!WriteBoundaryImages(scan, output_image_directory, "SA")) return 0;
	if (!WriteColorImages(scan, output_image_directory, "SA")) return 0;

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
			else if (!strcmp(*argv, "-overwrite")) { overwrite_detection = 1; }
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
		fprintf(stderr, "Usage: sadetect input_scene [output_directory] [-v] [-overwrite]\n");
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

















