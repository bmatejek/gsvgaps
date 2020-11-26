////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// Input/output
static char *input_run = NULL;
static int print_verbose = 0;
static int scanline = 0;
static double resolution = 1.0;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

static GSVScene *original_scene = NULL;
static GSVRun *original_run = NULL;
static GSVScene *aligned_scene = NULL;
static GSVRun *aligned_run = NULL;

// Useful constants
static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

static GSVScene *ReadScene(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate google scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Read scene 
  if (!scene->ReadFile(filename, true)) {
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
// Labeling functions
////////////////////////////////////////////////////////////////////////
static int autolabel(void) {
	// go through every image
    for (int is = 0; is < aligned_run->NSegments(); ++is) {
        // only do the first segment for now
        if (is != 0) continue;
        GSVSegment *aligned_segment = aligned_run->Segment(is);
        for (int ip = 0; ip < aligned_segment->NPanoramas(); ++ip) {
            if (ip < 146) continue;
            GSVPanorama *aligned_panorama = aligned_segment->Panorama(ip);
            for (int ii = 0; ii < aligned_panorama->NImages(); ++ii) {
                // only consider the three side images
                if (ii < 5 || ii == 8) continue;

                GSVImage *aligned_image = aligned_panorama->Image(ii);
                const R3Point& viewpoint = aligned_image->Pose().Viewpoint();

                // create output grid, initialize to 0
                R2Grid *detection_grid = new R2Grid(GSV_IMAGE_WIDTH / resolution, GSV_IMAGE_HEIGHT / resolution);
                R2Grid *label_grid = new R2Grid(GSV_IMAGE_WIDTH / resolution, GSV_IMAGE_HEIGHT / resolution);
                R2Grid *zbuffer_grid = new R2Grid(GSV_IMAGE_WIDTH / resolution, GSV_IMAGE_HEIGHT / resolution);
                for (int x = 0; x < detection_grid->XResolution(); ++x) {
                    for (int y = 0; y < detection_grid->YResolution(); ++y) {
                        detection_grid->SetGridValue(x, y, R2_GRID_UNKNOWN_VALUE);
                        label_grid->SetGridValue(x, y, R2_GRID_UNKNOWN_VALUE);
                        zbuffer_grid->SetGridValue(x, y, DBL_MAX);
                    }
                }

                // read in the super pixel image
                char image_name[4096];
                sprintf(image_name, "aligned/gsv_data/image_superpixels/%s/%02d_%06d_%02d_UndistortedImage.png", aligned_run->Name(), is, ip, ii);
                R2Grid superpixel_image;
                superpixel_image.Read(image_name);

                // go through each scanline
                for (int ia = 0; ia < aligned_segment->NScans(); ++ia) {
                    // for now skip these two scans
                    if (ia == 1) continue;
                    if (ia == 2) continue;

                    GSVScan *aligned_scan = aligned_segment->Scan(ia);
                    aligned_scan->ReadPoints();

                    R2Grid parse_detections, parse_labels;
                    R2Grid original_scanline, original_pointindex;

                    // get the parse detections
                    sprintf(image_name, "original/gsv_data/parse/%s/%02d_%02d_SA_ParseDetections.grd", original_run->Name(), is, ia);
                    parse_detections.Read(image_name);
                    sprintf(image_name, "original/gsv_data/parse/%s/%02d_%02d_SA_ParseLabels.grd", original_run->Name(), is, ia);
                    parse_labels.Read(image_name);

                    // get the scan and point indicies
                    sprintf(image_name, "original/gsv_data/laser_images/%s/%02d_%02d_SA_Scanline.grd", original_run->Name(), is, ia);
                    original_scanline.Read(image_name);
                    sprintf(image_name, "original/gsv_data/laser_images/%s/%02d_%02d_SA_PointIndex.grd", original_run->Name(), is, ia);
                    original_pointindex.Read(image_name);
                    // go through every point in the scan
                    for (int px = 0; px < original_scanline.XResolution(); ++px) {
                        for (int py = 0; py < original_scanline.YResolution(); ++py) {
                            if (parse_labels.GridValue(px, py) == R2_GRID_UNKNOWN_VALUE) continue;
                            RNScalar scanline_index_value = original_scanline.GridValue(px, py);
                            if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
                            RNScalar point_index_value = original_pointindex.GridValue(px, py);
                            if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;

                            int scanline_index = (int)(scanline_index_value + 0.5);
                            int point_index = (int)(point_index_value + 0.5);

                            GSVScanline *scanline = aligned_scan->Scanline(scanline_index);
                            const R3Point& position = scanline->PointPosition(point_index);

                            // project this point into the aligned image
                            const R2Point& middle = aligned_image->UndistortedPosition(position);

                            if (middle.X() < 0) continue;
                            if (middle.Y() < 0) continue;
                            if (middle.X() >= GSV_IMAGE_WIDTH) continue;
                            if (middle.Y() >= GSV_IMAGE_HEIGHT) continue;

                            // get distance between this point and viewpoint of camera to make sure it is close enough
                            RNScalar distance = R3Distance(position, viewpoint);
                            if (distance > zbuffer_grid->GridValue(middle.X() / resolution, middle.Y() / resolution)) continue;

                            // consider the four possible triangles
                            R3Point *bottom_left = NULL, *bottom_middle = NULL, *bottom_right = NULL;
                            R3Point *top_left = NULL, *top_middle = NULL, *top_right = NULL;
                           
                            // get the bottom left point if applicable
                            if (px != 0 && py != 0)  {
                                int pxo = px - 1;
                                int pyo = py - 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    bottom_left = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }

                            // get the bottom middle point if applicable
                            if (py != 0)  {
                                int pxo = px;
                                int pyo = py - 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    bottom_middle = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }

                            // get the bottom right point if applicable
                            if (px != original_scanline.XResolution() - 2 && py != 0)  {
                                int pxo = px + 1;
                                int pyo = py - 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    bottom_right = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }

                            // get the top left point if applicable
                            if (px != 0 && py != original_scanline.YResolution() - 2) {
                                int pxo = px - 1;
                                int pyo = py + 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    top_left = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }

                            // get the top middle point if applicable
                            if (py != original_scanline.YResolution() - 2) {
                                int pxo = px;
                                int pyo = py + 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    top_middle = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }
                            
                            // get the top right point if applicable
                            if (px != original_scanline.XResolution() - 2 && py != original_scanline.YResolution() - 2) {
                                int pxo = px + 1;
                                int pyo = py + 1;
                                RNScalar scanline_index_value_prime = original_scanline.GridValue(pxo, pyo);
                                RNScalar point_index_value_prime = original_pointindex.GridValue(pxo, pyo);
                                if (scanline_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    point_index_value_prime != R2_GRID_UNKNOWN_VALUE &&
                                    parse_detections.GridValue(px, py) == parse_detections.GridValue(pxo, pyo)) {
                                    int scanline_index_prime = (int)(scanline_index_value_prime + 0.5);
                                    int point_index_prime = (int)(point_index_value_prime + 0.5);
                                    GSVScanline *scanline_prime = aligned_scan->Scanline(scanline_index_prime);
                                    top_right = new R3Point(scanline_prime->PointPosition(point_index_prime));
                                }
                            }

                            R2Point *bottom_left_image = NULL, *bottom_middle_image = NULL, *bottom_right_image = NULL;
                            R2Point *top_left_image = NULL, *top_middle_image = NULL, *top_right_image = NULL;
                            // get the bottom left image point if applicable
                            if (bottom_left) {
                                const R2Point& point = aligned_image->UndistortedPosition(*bottom_left);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    bottom_left_image = new R2Point(point);
                                }
                            }

                            // get the bottom middle image point if applicable
                            if (bottom_middle) {
                                const R2Point& point = aligned_image->UndistortedPosition(*bottom_middle);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    bottom_middle_image = new R2Point(point);
                                }
                            }

                            // get the bottom right image point if applicable
                            if (bottom_right) {
                                const R2Point& point = aligned_image->UndistortedPosition(*bottom_right);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    bottom_right_image = new R2Point(point);
                                }
                            }

                            // get the top left image point if applicable
                            if (top_left) {
                                const R2Point& point = aligned_image->UndistortedPosition(*top_left);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    top_left_image = new R2Point(point);
                                }
                            }

                            // get the top middle image point if applicable
                            if (top_middle) {
                                const R2Point& point = aligned_image->UndistortedPosition(*top_middle);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    top_middle_image = new R2Point(point);
                                }
                            }
                            
                            // get the top right image point if applicable
                            if (top_right) {
                                const R2Point& point = aligned_image->UndistortedPosition(*top_right);
                                if (point.X() >= 0 && point.Y() >= 0 && point.X() < GSV_IMAGE_WIDTH && point.Y() < GSV_IMAGE_HEIGHT) {
                                    top_right_image = new R2Point(point);
                                }
                            }

                            int detection = (int)(parse_detections.GridValue(px, py) + 0.5);
                            int label = (int)(parse_labels.GridValue(px, py) + 0.5);
                            // consider the four sets of triangles
                            if (top_left_image && top_middle_image) {
                                detection_grid->RasterizeGridTriangleSetValue(*top_left_image / resolution, *top_middle_image / resolution, middle / resolution, detection);
                                label_grid->RasterizeGridTriangleSetValue(*top_left_image / resolution, *top_middle_image / resolution, middle / resolution, label);
                                zbuffer_grid->RasterizeGridTriangleSetValue(*top_left_image / resolution, *top_middle_image / resolution, middle / resolution, distance);
                            }
                            if (top_middle_image && top_right_image) {
                                detection_grid->RasterizeGridTriangleSetValue(*top_middle_image / resolution, *top_right_image / resolution, middle / resolution, detection);
                                label_grid->RasterizeGridTriangleSetValue(*top_middle_image / resolution, *top_right_image / resolution, middle / resolution, label);
                                zbuffer_grid->RasterizeGridTriangleSetValue(*top_middle_image / resolution, *top_right_image / resolution, middle / resolution, distance);
                            }
                            if (bottom_left_image && bottom_middle_image) {
                                detection_grid->RasterizeGridTriangleSetValue(*bottom_left_image / resolution, *bottom_middle_image / resolution, middle / resolution, detection);
                                label_grid->RasterizeGridTriangleSetValue(*bottom_left_image / resolution, *bottom_middle_image / resolution, middle / resolution, label);
                                zbuffer_grid->RasterizeGridTriangleSetValue(*bottom_left_image / resolution, *bottom_middle_image / resolution, middle / resolution, distance);
                            }
                            if (bottom_middle_image && bottom_right_image) {
                                detection_grid->RasterizeGridTriangleSetValue(*bottom_middle_image / resolution, *bottom_right_image / resolution, middle / resolution, detection);
                                label_grid->RasterizeGridTriangleSetValue(*bottom_middle_image / resolution, *bottom_right_image / resolution, middle / resolution, label);                            
                                zbuffer_grid->RasterizeGridTriangleSetValue(*bottom_middle_image / resolution, *bottom_right_image / resolution, middle / resolution, distance);
                            }
                            // delete all the points
                            delete bottom_right; delete bottom_middle; delete bottom_left;
                            delete top_right; delete top_middle; delete top_left;
                            delete bottom_left_image; delete bottom_middle_image; delete bottom_right_image;
                            delete top_left_image; delete top_middle_image; delete top_right_image;
                        }
                    }
                    aligned_scan->ReleasePoints();
                }
                // determine the number of super pixels in the image
                int numsuperpixels = superpixel_image.Maximum() + 1;
                int numlabels = label_grid->Maximum() + 2;
                if (numlabels < 0) numlabels = 1;
                int numdetections = detection_grid->Maximum() + 2;
                if (numdetections < 0) numdetections = 1;
                vector<vector<int> > pixel_labels = vector<vector<int> >(numsuperpixels);
                vector<vector<int> > pixel_detections = vector<vector<int> >(numsuperpixels);
                vector<int> superpixel_size = vector<int>(numsuperpixels);
                
                for (int nsp = 0; nsp < numsuperpixels; ++nsp) {
                    pixel_labels[nsp] = vector<int>(numlabels);
                    pixel_detections[nsp] = vector<int>(numdetections);
                    for (int nl = 0; nl < numlabels; ++nl) {
                        pixel_labels[nsp][nl] = 0;
                    }
                    for (int nd = 0; nd < numdetections; ++nd) {
                        pixel_detections[nsp][nd] = 0;
                    }
                }

                // go through the superpixel image and determine the best labels and detections
                for (int x = 0; x < GSV_IMAGE_WIDTH / resolution; ++x) {
                    for (int y = 0; y < GSV_IMAGE_HEIGHT / resolution; ++y) {
                        int superpixel = superpixel_image.GridValue(x, y);
                        int detection = detection_grid->GridValue(x, y);
                        int label = label_grid->GridValue(x, y);
                        superpixel_size[superpixel]++;
                        if (detection == R2_GRID_UNKNOWN_VALUE)
                            pixel_detections[superpixel][numdetections - 1]++;
                        else
                            pixel_detections[superpixel][detection]++;

                        if (label == R2_GRID_UNKNOWN_VALUE)
                            pixel_labels[superpixel][numlabels - 1]++;
                        else
                            pixel_labels[superpixel][label]++;

                    }
                }

                R2Grid output_detection = R2Grid(GSV_IMAGE_WIDTH / resolution, GSV_IMAGE_HEIGHT / resolution);
                R2Grid output_label = R2Grid(GSV_IMAGE_WIDTH / resolution, GSV_IMAGE_HEIGHT / resolution);

                // each superpixel gets the largest vote getter
                vector<int> superpixel_label_map = vector<int>(numsuperpixels);
                vector<int> superpixel_detection_map = vector<int>(numsuperpixels);
                
                for (int nsp = 0; nsp < numsuperpixels; ++nsp) {
                    superpixel_label_map[nsp] = numlabels - 1;
                }
                for (int nsp = 0; nsp < numsuperpixels; ++nsp) {
                    superpixel_detection_map[nsp] = numdetections - 1;
                }

                for (int nsp = 0; nsp < numsuperpixels; ++nsp) {
                    int mostlabeled = 0;
                    int mostdetected = 0;
                    for (int nl = 0; nl < numlabels; ++nl) {
                        if (pixel_labels[nsp][nl] > mostlabeled && !(nl == numlabels - 1 && (pixel_labels[nsp][nl] / (double)mostlabeled < 3.0))) {
                            mostlabeled = pixel_labels[nsp][nl];
                            if (pixel_labels[nsp][nl] / (double)superpixel_size[nsp] < 0.05)
                                superpixel_label_map[nsp] = numlabels - 1;
                            else
                                superpixel_label_map[nsp] = nl;
                        }
                    }
                    for (int nd = 0; nd < numdetections; ++nd) {
                        if (pixel_detections[nsp][nd] > mostdetected && !(nd == numdetections - 1 && (pixel_detections[nsp][nd] / (double)mostdetected < 3.0))) {
                            mostdetected = pixel_detections[nsp][nd];
                            if (pixel_detections[nsp][nd] / (double)superpixel_size[nsp] < 0.05)
                                superpixel_detection_map[nsp] = numdetections - 1;
                            else
                                superpixel_detection_map[nsp] = nd;
                        }
                    }
                }

                for (int x = 0; x < GSV_IMAGE_WIDTH / resolution; ++x) {
                    for (int y = 0; y < GSV_IMAGE_HEIGHT / resolution; ++y) {
                        int superpixel = superpixel_image.GridValue(x, y);
                        if (superpixel_label_map[superpixel] == numlabels - 1)
                            output_label.SetGridValue(x, y, R2_GRID_UNKNOWN_VALUE);
                        else
                            output_label.SetGridValue(x, y, superpixel_label_map[superpixel]);
                        if (superpixel_detection_map[superpixel] == numdetections - 1)
                            output_detection.SetGridValue(x, y, R2_GRID_UNKNOWN_VALUE);
                        else
                            output_detection.SetGridValue(x, y, superpixel_detection_map[superpixel]);
                    }
                }

                printf("%02d_%06d_%02d\n", is, ip, ii);
                // output the original image files
                sprintf(image_name, "aligned/gsv_data/labeled_images/%s/%02d_%06d_%02d_Detections.grd", aligned_run->Name(), is, ip, ii);
                output_detection.Write(image_name);
                sprintf(image_name, "aligned/gsv_data/labeled_images/%s/%02d_%06d_%02d_Labels.grd", aligned_run->Name(), is, ip, ii);
                output_label.Write(image_name);

                // output the original files
                /*sprintf(image_name, "aligned/gsv_data/labeled_images/%s/%02d_%06d_%02d_DetectionsX.grd", aligned_run->Name(), is, ip, ii);
                detection_grid->Write(image_name);
                sprintf(image_name, "aligned/gsv_data/labeled_images/%s/%02d_%06d_%02d_LabelsX.grd", aligned_run->Name(), is, ip, ii);
                label_grid->Write(image_name);
                sprintf(image_name, "aligned/gsv_data/labeled_images/%s/%02d_%06d_%02d_Distances.grd", aligned_run->Name(), is, ip, ii);
                zbuffer_grid->Write(image_name);
                */
                delete detection_grid;
                delete label_grid;
                delete zbuffer_grid;
            }
        }
    }

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage(void) {
	fprintf(stderr, "Usage: autolabel run [-v]\n");

	// return failure status
	return 0;
}

static int ParseArgs(int argc, char **argv) {
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
	  else if (!strcmp(*argv, "-resolution")) { argv++; argc--; resolution = atof(*argv); }
	  else if (!strcmp(*argv, "-scanline")) { argv++; argc--; scanline = atoi(*argv); }
	  else if (!strcmp(*argv, "-usage")) { return Usage(); }
      else { 
        fprintf(stderr, "Invalid program argument: %s\n", *argv); 
		return Usage();
      }
      argv++; argc--;
    }
    else {
      if (!input_run) input_run = *argv;
	  else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_run) return Usage();

  // Return OK status 
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	RNTime start_time;
	start_time.Read();

	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

    char input_scene_name[4096];
	// Read the original scene
	sprintf(input_scene_name, "original/gsv_data/laser_points/%s/%s.gsv", input_run, input_run);
	original_scene = ReadScene(input_scene_name);
	if (!original_scene) exit(-1);
    original_run = original_scene->Run(0);

    // Read the aligned scene
    sprintf(input_scene_name, "aligned/gsv_data/laser_points/%s/%s.gsv", input_run, input_run);
    aligned_scene = ReadScene(input_scene_name);
    if (!aligned_scene) exit(-1);
    aligned_run = aligned_scene->Run(0);

	// label all of the images
	autolabel();

	if (print_verbose) {
		printf("Populated auto labels ... \n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
	}

	// Return success 
	return 0;
}


