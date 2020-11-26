// Source file for the google viewer program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static const char *input_scene_name = NULL;
static const char *output_image_directory = NULL;
static int capture_DA_images = 1;
static int print_verbose = 0;

static RNLength DA_image_spacing = 0.1;

////////////////////////////////////////////////////////////////////////
// Grid names
////////////////////////////////////////////////////////////////////////

// Grid names
const char *base_image_names[] = {
    "ParseLabels"
};

int num_base_image_names = sizeof(base_image_names) / sizeof(const char *);

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
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
// Utility types and functions
////////////////////////////////////////////////////////////////////////

static void
InterpolateMissingColumns(R2Grid& grid)
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
                    if (abs(ix0 - ix1) > 20) continue;
                    RNScalar t = (double)(ix - ix0) / (double)(ix1 - ix0);
                    RNScalar value = value0; // (1 - t)*value0 + t*value1;
                    grid.SetGridValue(ix, iy, value);
                }
            }
        }

        // Remember last column with values
        ix0 = ix1;
    }
}


////////////////////////////////////////////////////////////////////////
// DA images (Distance vs. Angle)
////////////////////////////////////////////////////////////////////////

static int
WriteDAImage(GSVScan *scan, const char *output_image_directory,
const R2Grid& sa_travel_distance_grid, const R2Grid& sa_viewpoint_distance_grid,
const char *grid_name)
{
    // Get convenient variables
    if (!scan) return 0;
    if (scan->NScanlines() == 0) return 0;
    int scan_index = scan->SegmentIndex();
    GSVSegment *segment = scan->Segment();
    if (!segment) return 0;
    int segment_index = segment->RunIndex();
    GSVRun *run = segment->Run();
    if (!run) return 0;

    // Skip if already done
    char da_name[4096];
    sprintf(da_name, "%s/%s/%02d_%02d_DA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
    if (RNFileExists(da_name)) return 1;

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



static int
WriteDAImages(GSVScan *scan, const char *output_image_directory)
{
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

static int WriteImages(GSVScene *scene, const char *output_image_directory)
{
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
            //if (is != 0) continue;
            for (int ia = 0; ia < segment->NScans(); ia++) {
                if (ia == 1) continue;
                GSVScan *scan = segment->Scan(ia);
                if (capture_DA_images && !WriteDAImages(scan, output_image_directory)) return 0;
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
        fprintf(stderr, "Usage: sa2da input_scene [output_directory] [options]\n");
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

int main(int argc, char **argv)
{
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

















