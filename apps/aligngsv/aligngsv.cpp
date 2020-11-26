// Source file for the google street view pose optimization program 



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "R3CatmullRomSpline.h"
#include "fglut/fglut.h"
#include "poseoptimization.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// Program arguments

static char *input_scene_name = NULL;
static char *output_scene_name = NULL;
static char *input_correspondences_name = NULL;
static char *output_correspondences_name = NULL;
static char *input_transformations_name = NULL;
static char *output_transformations_name = NULL;
static char *output_error_plot_name = NULL;
static int run_count = 0;
static int segment_count = 0;
static int create_scan_pole_features = 0;
static int create_scan_curb_features = 0;
static int create_scan_edge_features = 0;
static int create_scan_plane_features = 0;
static int create_image_corner_features = 0;
static int create_image_sift_features = 0;
static int create_image_line_features = 0;
static int create_features = 0;
static int create_point_point_correspondences = 0;
static int create_plane_plane_correspondences = 0;
static int create_pixel_point_correspondences = 0;
static int create_pixel_plane_correspondences = 0;
static int create_pixel_pixel_correspondences = 0;
static int create_correspondences = 0;
static int create_transformations = 0;
static int update_pixel_features = 0;
static int load_points_from_images = 0;
static int load_points_dynamically = 1;
static int print_verbose = 0;
static int print_debug = 0;
static int interactive = 0;



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Data variables

static GSVScene *scene = NULL;
static GSVPoseOptimization *optimization = NULL;
static GSVFeatureCorrespondence *last_correspondence = NULL;
static GSVFeature *last_feature = NULL;
static int num_initial_correspondences = 0;



// Display variables

static int show_scans[2] = { 1, 1 };
static int show_images[2] = { 0, 0 };
static int show_meshes[2] = { 0, 0 };
static int show_cameras[2] = { 0, 0 };
static int show_lasers[2] = { 1, 1 };
static int show_paths[2] = { 0, 0 };
static int show_bboxes[2] = { 0, 0 };
static int show_features[2] = { 0, 0 };
static int show_segments[2] = { 0, 0 };
static int show_correspondences[2] = { 1, 1 };
static int show_optimization[2] = { 1, 1 };
static R3Viewer *viewer[2] = { NULL, NULL };
static GSVImage *selected_image[2] = { NULL, NULL };
static GSVScanline *selected_scanline[2] = { NULL, NULL };
static GSVFeature *selected_feature[2] = { NULL, NULL };
static R3Point center_point[2] = { R3Point(0,0,0), R3Point(0,0,0) };
static int viewpoint_scheme[2] = { 0, 0 };
static int color_scheme[2] = { 1, 1 };
static float point_size[2] = { 1, 1 };
static int timestamp_radius = 10;
static double image_plane_distance = 20;
static double rigidity = 1.0;
static int solver = RN_POLYNOMIAL_CERES_SOLVER;
static int split_screen = 0;

// limit the number of runs shown
static int show_all_runs = 1;
static int run_number = 0;

// limit the number of segments shown
static int show_all_segments = 1;
static int segment_number = 0;

// total number of intersections found
static int total_intersections = 0;


// Useful constants

static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 1024;
static int GLUTwindow_width = 2 * GLUTwindow_height * GSV_IMAGE_WIDTH / GSV_IMAGE_HEIGHT;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;
static int GLUTside = 0;



// Color schemes

enum {
  NO_COLOR_SCHEME,
  SEGMENT_COLOR_SCHEME,
  VIEWPOINT_DEPTH_COLOR_SCHEME,
  HEIGHT_COLOR_SCHEME,
  PICK_COLOR_SCHEME
};



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
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
  if (!scene->ReadFile(filename, !load_points_dynamically)) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    
    int totalRuns = 0;
    int totalCameras = 0;
    int totalLasers = 0;
    int totalSegments = 0;
    int totalScans = 0;
    int totalTapestries = 0;
    int totalPanoramas = 0;
    int totalImages = 0;
    int totalScanlines = 0;
    
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      if (show_all_runs || ir == run_number) {
       GSVRun *run = scene->Run(ir);
       bool containsSegment = false;
       printf("Total number of segments: %d\n", run->NSegments());
       for (int is = 0; is < run->NSegments(); is++) {
         if (show_all_segments || is == segment_number) {
           GSVSegment *segment = run->Segment(is);
           printf("Segment Number %d - Number of Images %d\n", is, segment->NImages());
           containsSegment = true;
           totalSegments++;
           totalScans += segment->NScans();
           totalTapestries += segment->NTapestries();
           totalPanoramas += segment->NPanoramas();
           totalImages += segment->NImages();
           totalScanlines += segment->NScanlines();
         }
       }
       if (containsSegment) {
         totalRuns++;
         totalCameras += run->NCameras();
         totalLasers += run->NLasers();
       }
     }
   }

   printf("  # Runs = %d\n", totalRuns);
   printf("  # Cameras = %d\n", totalCameras);
   printf("  # Lasers = %d\n", totalLasers);
   printf("  # Segments = %d\n", totalSegments);
   printf("  # Scans = %d\n", totalScans);
   printf("  # Tapestries = %d\n", totalTapestries);
   printf("  # Panoramas = %d\n", totalPanoramas);
   printf("  # Images = %d\n", totalImages);
   printf("  # Scanlines = %d\n", totalScanlines);
   fflush(stdout);

    /* code before allowing selections of runs and segments */
   if (0) {
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
}

  // Return scene
return scene;
}



static int
WriteScene(GSVScene *scene, GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Apply pose transformations
  if (optimization) {
    if (!optimization->ApplyOptimizedTransformationsToScene()) return 0;
  }

  // Write scene
  if (!scene->WriteFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote scene to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    
    int totalRuns = 0;
    int totalCameras = 0;
    int totalLasers = 0;
    int totalSegments = 0;
    int totalScans = 0;
    int totalTapestries = 0;
    int totalPanoramas = 0;
    int totalImages = 0;
    int totalScanlines = 0;

    for (int ir = 0; ir < scene->NRuns(); ir++) {
      if (show_all_runs || ir == run_number) {
       GSVRun *run = scene->Run(ir);
       bool containsSegment = false;
       for (int is = 0; is < run->NSegments(); is++) {
         if (show_all_segments || is == segment_number) {
           GSVSegment *segment = run->Segment(is);
           containsSegment = true;
           totalSegments++;
           totalScans += segment->NScans();
           totalTapestries += segment->NTapestries();
           totalPanoramas += segment->NPanoramas();
           totalImages += segment->NImages();
           totalScanlines += segment->NScanlines();
         }
       }
       if (containsSegment) {
         totalRuns++;
         totalCameras += run->NCameras();
         totalLasers += run->NLasers();	
       }
     }
   }

   printf("  # Runs = %d\n", totalRuns);
   printf("  # Cameras = %d\n", totalCameras);
   printf("  # Lasers = %d\n", totalLasers);
   printf("  # Segments = %d\n", totalSegments);
   printf("  # Scans = %d\n", totalScans);
   printf("  # Tapestries = %d\n", totalTapestries);
   printf("  # Panoramas = %d\n", totalPanoramas);
   printf("  # Images = %d\n", totalImages);
   printf("  # Scanlines = %d\n", totalScanlines);
   fflush(stdout);

    /* code before allowing selections of runs and segments */
   if (0) {
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
}

  // Return success
return 1;
}



static int
ReadCorrespondences(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read correspondences
  if (!optimization->ReadFile(filename)) return 0;

  // Remember number of initial correspondences
  num_initial_correspondences = optimization->NCorrespondences();

  // Print statistics
  if (print_verbose) {
    printf("Read correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  Score = %g\n", optimization->Score());
    printf("  RMSD = %g\n", optimization->RMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteCorrespondences(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  
  // Determine if should apply transformations when write
  RNBoolean apply_pose_transformations = (output_scene_name) ? TRUE : FALSE;

  // Write correspondences
  if (!optimization->WriteFile(filename, apply_pose_transformations)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote correspondences to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  Score = %g\n", optimization->Score());
    printf("  RMSD = %g\n", optimization->RMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadTransformations(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read transformations
  if (!optimization->ReadTransformationsFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Read transformations from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Score = %g\n", optimization->Score());
    printf("  RMSD = %g\n", optimization->RMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteTransformations(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write transformations
  if (!optimization->WriteTransformationsFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote transformations to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Score = %g\n", optimization->Score());
    printf("  RMSD = %g\n", optimization->RMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static R2Image *
ReadUndistortedImage(GSVImage *image)
{
  // Get convenient variables
  int image_index = image->PanoramaIndex();
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return NULL;
  int panorama_index = panorama->RunIndex();
  GSVSegment *segment = panorama->Segment();
  if (!segment) return NULL;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return NULL;
  GSVScene *scene = run->Scene();
  if (!scene) return NULL;
  R2Image *img = NULL;

  // Construct scaled image name
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%06d_%02d_UndistortedImage_2X.jpg", 
    scene->CacheDataDirectoryName(), run->Name(), segment_index, panorama_index, image_index);
  
  // Check if scaled image file exists
  if (RNFileExists(image_name)) {
    // Allocate image
    img = new R2Image();
    if (!img) {
      fprintf(stderr, "Unable to allocate image for %s\n", image_name);
      return NULL;
    }

    // Read file
    if (!img->Read(image_name)) {
      fprintf(stderr, "Unable to read image from %s\n", image_name);
      delete img;
      return NULL;
    }
  }
  else {
    // Read full size image
    img = image->UndistortedImage();
    if (!img) return NULL;
  }

  // Return undistorted image
  return img;
}



static int
WriteErrorPlot(GSVPoseOptimization *optimization, const char *filename = NULL)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Create histogram
  double bin_spacing = 0.25;
  RNLength max_error = 10;
  int nbins = (int) (max_error / bin_spacing);
  double *bins = new double [nbins];
  for (int i = 0; i < nbins; i++) bins[i] = 0;
    for (int i = 0; i < optimization->NCorrespondences(); i++) {
      GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
      GSVFeature *feature0 = correspondence->Feature(0);
      GSVFeature *feature1 = correspondence->Feature(1);
      if (feature0->feature_type != GSV_SCAN_POINT_FEATURE_TYPE) continue;
      if (feature1->feature_type != GSV_SCAN_POINT_FEATURE_TYPE) continue;
      R3Point position0 = optimization->FeaturePosition(feature0, TRUE);
      R3Point position1 = optimization->FeaturePosition(feature1, TRUE);
      RNLength error = R3Distance(position0, position1);
      double bin = error / bin_spacing;
      int bin0 = (int) bin;
      int bin1 = bin0 + 1;
      if (bin0 >= nbins) bin0 = nbins - 1;
      if (bin1 >= nbins) bin1 = nbins - 1;
      double t = bin - bin0;
      bins[bin0] += 1-t;
      bins[bin1] += t;
    }

  // Open file
    FILE *fp = stdout;
    if (filename) fp = fopen(filename, "w");
    if (!fp) { 
      fprintf(stderr, "Unable to open error plot %s\n", filename); 
      return 0; 
    }

  // Write CDF error plot
    double sum = 0;
    double total = optimization->NCorrespondences();
    if (total == 0) total = 1;
    for (int i = 0; i < nbins; i++) {
      sum += bins[i];
      fprintf(fp, "%g %g\n", (i+0.5) * bin_spacing, sum / total);
    }

  // Close file
    if (filename) fclose(fp);

  // Delete histogram 
    delete [] bins;

  // Print statistics
    if (print_verbose) {
      printf("Wrote error plot to %s ...\n", filename);
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Features = %d\n", optimization->NFeatures());
      printf("  # Correspondences = %d\n", optimization->NCorrespondences());
      printf("  Score = %g\n", optimization->Score());
      printf("  RMSD = %g\n", optimization->RMSD());
      fflush(stdout);
    }

  // Return success
    return 1;
  }



////////////////////////////////////////////////////////////////////////
// Feature creation functions
////////////////////////////////////////////////////////////////////////

  static int
  CreateFeatures(GSVPoseOptimization *optimization)
  {
  // Start statistics
    RNTime start_time;
    start_time.Read();
    int saved_nfeatures = optimization->NFeatures();
    if (print_verbose) {
      printf("Creating features ...\n");
      fflush(stdout);
    }

  // Create scan pole features
    if (create_scan_pole_features) {
      if (!optimization->CreateScanPoleFeatures()) return 0;
    }
    
  // Create scan curb features
    if (create_scan_curb_features) {
      if (!optimization->CreateScanCurbFeatures()) return 0;
    }
    
  // Create scan edge features
    if (create_scan_edge_features) {
      if (!optimization->CreateScanEdgeFeatures()) return 0;
    }
    
  // Create scan plane features
    if (create_scan_plane_features) {
      if (!optimization->CreateScanPlaneFeatures()) return 0;
    }
    
  // Create image corner features
    if (create_image_corner_features) {
      if (!optimization->CreateImageCornerFeatures()) return 0;
    }

  // Create image sift features
    if (create_image_sift_features) {
      printf("Here\n");
      if (!optimization->CreateImageSiftFeatures()) return 0;
    }

  // Create image line features
    if (create_image_line_features) {
      if (!optimization->CreateImageLineFeatures()) return 0;
    }

  // Print statistics
    if (print_verbose) {
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # New Features = %d\n", optimization->NFeatures() - saved_nfeatures);
      fflush(stdout);
    }

  // Return success
    return 1;
  }



////////////////////////////////////////////////////////////////////////
// Temporary
////////////////////////////////////////////////////////////////////////

  static int 
  UpdatePixelFeatures(GSVPoseOptimization *optimization, const R3Box& bbox)
  {
  // Start statistics
    RNTime start_time;
    start_time.Read();
    if (print_verbose) {
      printf("Updating pixel features ...\n");
      fflush(stdout);
    }

  // Delete features that are not in bbox
    int nfeatures1 = optimization->NFeatures();
    if (!bbox.IsEmpty()) {
      RNArray<GSVFeature *> features = optimization->features.features;
      for (int i = 0; i < features.NEntries(); i++) {
        GSVFeature *feature = features.Kth(i);
        if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
        const R3Point& viewpoint = feature->image->Pose().Viewpoint();
        if (!R3Contains(bbox, viewpoint)) {
          optimization->RemoveFeature(feature); 
        }
      }
    }

  // Delete features that are not in center of image
    int nfeatures2 = optimization->NFeatures();
    RNScalar center_size = 0.8;
    if (center_size < 1.0) {
      RNScalar lo = 0.5 * (1 - center_size);
      RNScalar hi = 1 - lo;
      RNArray<GSVFeature *> features = optimization->features.features;
      for (int i = 0; i < features.NEntries(); i++) {
        GSVFeature *feature = features.Kth(i);
        if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
        RNScalar x = feature->image_position.X() / feature->image->Width();
        RNScalar y = feature->image_position.X() / feature->image->Width();
        if ((x < lo) || (x > hi) || (y < lo) || (y > hi)) {
          optimization->RemoveFeature(feature); 
        }
      }
    }

  // Update positions of pixel features
    if (!optimization->UpdatePixelPositionsFromRayIntersections()) exit(-1);

  // Delete features with no mesh intersection
    int nfeatures3 = optimization->NFeatures();
    if (1) {
      RNArray<GSVFeature *> features = optimization->features.features;
      for (int i = 0; i < features.NEntries(); i++) {
        GSVFeature *feature = features.Kth(i);
        if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
        if (feature->image_t > 0.0) continue;
        optimization->RemoveFeature(feature);
      }
    }

  // Print statistics
    if (print_verbose) {
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Features = %d\n", optimization->NFeatures());
      printf("  # Features outside box = %d\n", nfeatures1 - nfeatures2);
      printf("  # Features near border = %d\n", nfeatures2 - nfeatures3);
      printf("  # Features not on mesh = %d\n", nfeatures3 - optimization->NFeatures());
      fflush(stdout);
    }

  // Return success
    return 1;
  }



////////////////////////////////////////////////////////////////////////
// Correspondence creation functions
////////////////////////////////////////////////////////////////////////

  static int
  CreateCorrespondences(GSVPoseOptimization *optimization)
  {
  // Start statistics
    RNTime start_time;
    start_time.Read();
    int saved_ncorrespondences = optimization->NCorrespondences();
    if (print_verbose) {
      printf("Creating correspondences ...\n");
      fflush(stdout);
    }

  // Reset correspondences to initial set
    optimization->correspondences.correspondences.Truncate(num_initial_correspondences);

  // Select which correspondences to create
    optimization->create_point_point_correspondences = create_point_point_correspondences;
    optimization->create_plane_plane_correspondences = create_plane_plane_correspondences;
    optimization->create_pixel_point_correspondences = create_pixel_point_correspondences;
    optimization->create_pixel_plane_correspondences = create_pixel_plane_correspondences;
    optimization->create_pixel_pixel_correspondences = create_pixel_pixel_correspondences;

  // Update optimization
    if (!optimization->CreateCorrespondences()) return 0;

  // Print statistics
    if (print_verbose) {
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # New Correspondences = %d\n", optimization->NCorrespondences() - saved_ncorrespondences);
      printf("  Score = %g\n", optimization->Score());
      printf("  RMSD = %g\n", optimization->RMSD());
      fflush(stdout);
    }

  // Return success
    return 1;
  }



////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Optimization functions
////////////////////////////////////////////////////////////////////////

  static int
  Solve(GSVPoseOptimization *optimization, 
    RNBoolean update_translations = TRUE, RNBoolean update_rotations = TRUE, 
    RNBoolean update_laser_transformations = TRUE, RNBoolean update_camera_transformations = TRUE, 
    RNBoolean update_path_transformations = TRUE, RNBoolean update_pixel_depths = TRUE,
    RNBoolean include_scan_scan_correspondences = TRUE, 
    RNBoolean include_image_image_correspondences = TRUE, 
    RNBoolean include_scan_image_correspondences = TRUE,
    RNScalar rigidity = 0.5, int solver = 0)
  {
  // Start statistics
    RNTime start_time;
    start_time.Read();

  // Update optimization
    if (!optimization->Solve(update_translations, update_rotations,
      update_laser_transformations, update_camera_transformations, 
      update_path_transformations, update_pixel_depths,
      include_scan_scan_correspondences, include_image_image_correspondences, include_scan_image_correspondences,
      rigidity, solver)) {
      fprintf(stderr, "Unable to update optimization\n");
    return 0;
  }
  
  // Print statistics
  if (print_verbose) {
    printf("Optimized pose transformations ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  Score = %g\n", optimization->Score());
    printf("  RMSD = %g\n", optimization->RMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
RandomizePoseTranslations(GSVPoseOptimization *optimization)
{
  // Assign random translation to every path
  // for (int i = 0; i < optimization->path_vertices.NEntries(); i++) {
    // GSVPathVertex *vertex = optimization->path_vertices.Kth(i);
    // vertex->rotation = R3zero_vector;
    // vertex->translation.Reset(50.0 * RNRandomScalar(), 50.0 * RNRandomScalar(), 0);
  // }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Scanline data functions
////////////////////////////////////////////////////////////////////////

static int
ReadScanData(GSVScan *scan, RNScalar timestamp) 
{
  // Check scan 
  if (!scan) return 1;
  ScanData *scan_data = (ScanData *) scan->Data();
  if (scan_data) { scan_data->timestamp = timestamp; return 1; }
  if (scan->NScanlines() == 0) return 1;

  // Get convenient variables
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  GSVRun *run = segment->Run();
  if (!run) return 0;
  int segment_index = segment->RunIndex();
  char image_name[4096];

  // Read points
  if (!scan->ReadPoints()) {
    fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Fill scan data
  scan_data = new ScanData();
  scan->SetData(scan_data);
  scan_data->scan = scan;
  scan_data->mesh = NULL; 
  // scan_data->mesh = scan->Mesh();
  scan_data->timestamp = timestamp;
  scan_data->index = scan->SceneIndex();
  if (load_points_from_images) {
    static const char *image_directory = "gsv_data/laser_images";
    sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_scanline_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_point_index_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_position_x_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_position_y_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_position_z_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentId.grd", image_directory, run->Name(), segment_index, scan_index);
    if (!scan_data->da_plane_segment_id_grid.Read(image_name)) return 0;
  }

  // Return success
  return 1;
}



static int
ReleaseScanData(GSVScan *scan)
{
  // Check scan
  if (!scan) return 1;
  if (!scan->Data()) return 1;

  // Release points
  if (!scan->ReleasePoints()) {
    GSVSegment *segment = scan->Segment();
    GSVRun *run = segment->Run();
    fprintf(stderr, "Unable to release points for %s %02d %02d\n", run->Name(), segment->RunIndex(), scan->SegmentIndex());
    return 0;
  }

  // Delete scan data
  ScanData *scan_data = (ScanData *) scan->Data();
  if (scan_data->mesh) delete scan_data->mesh;
  delete scan_data;
  scan->SetData(NULL);

  // Return success
  return 1;
}



static int
ToggleScanData(GSVSegment *segment, RNScalar timestamp) 
{
  // Check segment
  if (!segment) return 1;

  // Fill scan data
  for (int i = 0; i < segment->NScans(); i++) {
    GSVScan *scan = segment->Scan(i);
    if (scan->Data()) ReleaseScanData(scan);
    else ReadScanData(scan, timestamp);
  }      

  // Return success
  return 1;
}



static int
SelectScanData(GSVSegment *selected_segment, RNScalar timestamp) 
{
  // Fill scan data for segment, and empty all others
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           if (segment == selected_segment) ReadScanData(scan, timestamp);
           else ReleaseScanData(scan);
         }
       }
     }
   }
 }

  // Return success
 return 1;
}



////////////////////////////////////////////////////////////////////////
// Transformation functions
////////////////////////////////////////////////////////////////////////

static GSVScanline *pushed_scanline = NULL;
static GSVImage *pushed_image = NULL;



static void
PushTransformation(GSVScanline *scanline)
{
  // Check pushed scanline
  if (pushed_scanline) return;
  if (pushed_image) return;

  // Set transformation
  optimization->OptimizedTransformation(scanline).Push();

  // Remember pushed scanline
  pushed_scanline = scanline;
  pushed_image = NULL;
}



static void
PushTransformation(GSVImage *image)
{
  // Check pushed scanline
  if (pushed_scanline) return;
  if (pushed_image) return;
  
  // Set transformation
  optimization->OptimizedTransformation(image).Push();
  
  // Remember pushed scanline
  pushed_image = image;
  pushed_scanline = NULL;
}




static void
PopTransformation(void)
{
  // Check pushed scanline
  if (!pushed_scanline && !pushed_image) return;

  // Pop transformation 
  // if (pushed_scanline) optimization->OptimizedTransformation(pushed_scanline).Pop();
  // else if (pushed_image) optimization->OptimizedTransformation(pushed_image).Pop();
  R3identity_affine.Pop(); // (a bit of a hack for efficiency)

  // Update pushed scanline
  pushed_scanline = NULL;
  pushed_image = NULL;
}



////////////////////////////////////////////////////////////////////////
// Color functions
////////////////////////////////////////////////////////////////////////

static void
LoadSegmentColor(GSVSegment *segment, int color_scheme)
{
  // Check segment
  if (!segment) return;

  // Segment colors
  const int max_segment_colors = 12;
  const RNRgb segment_colors[max_segment_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),

  };

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    // Find run index and segment index
    GSVRun *run = segment->Run();
    if (!run) glColor3d(0.0, 0.0, 0.0);
    int segment_index = segment->RunIndex();
    int run_index = run->SceneIndex();

    // Load segment color
    RNLoadRgb(segment_colors[(7*run_index + segment_index) % max_segment_colors]);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int segment_index = segment->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (segment_index >> 16) & 0xFF;
    rgba[1] = (segment_index >> 8) & 0xFF;
    rgba[2] = segment_index & 0xFF;
    rgba[3] = 0xFE;
    glColor4ubv(rgba);
  }
}



static void
LoadScanlineColor(GSVScanline *scanline, int side, int color_scheme)
{
  // Check scanline
  if (!scanline) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    GSVScan *scan = scanline->Scan();
    if (!scan) glColor3d(0.0, 0.0, 0.0);
    GSVSegment *segment = scan->Segment();
    LoadSegmentColor(segment, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int scanline_index = 0;
    ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
    if (scanline_data) scanline_index = scanline_data->index;
    else scanline_index = scanline->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (scanline_index >> 16) & 0xFF;
    rgba[1] = (scanline_index >> 8) & 0xFF;
    rgba[2] = scanline_index & 0xFF;
    rgba[3] = 0xFD;
    glColor4ubv(rgba);
  }
}



static void
LoadImageColor(GSVImage *image, int side, int color_scheme)
{
  // Check image
  if (!image) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    if (image == selected_image[side]) { glColor3d(1.0, 1.0, 1.0); return; }
    GSVPanorama *panorama = image->Panorama();
    if (!panorama) { glColor3d(0.0, 0.0, 0.0); return; }
    GSVSegment *segment = panorama->Segment();
    LoadSegmentColor(segment, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int image_index = 0;
    ImageData *image_data = (ImageData *) image->Data();
    if (image_data) image_index = image_data->index;
    else image_index = image->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (image_index >> 16) & 0xFF;
    rgba[1] = (image_index >> 8) & 0xFF;
    rgba[2] = image_index & 0xFF;
    rgba[3] = 0xFC;
    glColor4ubv(rgba);
  }
}



static void
LoadFeatureColor(GSVFeature *feature, int side, int color_scheme, int feature_index)
{
  // Check feature
  if (!feature) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    if (feature == selected_feature[side]) { glColor3d(1.0, 1.0, 1.0); return; }
    GSVScanline *scanline = feature->scanline;
    if (!scanline) { glColor3d(0.0, 0.0, 0.0); return; }
    LoadScanlineColor(scanline, side, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    unsigned char rgba[4];
    rgba[0] = (feature_index >> 16) & 0xFF;
    rgba[1] = (feature_index >> 8) & 0xFF;
    rgba[2] = feature_index & 0xFF;
    rgba[3] = 0xFB;
    glColor4ubv(rgba);
  }
}



static void
LoadHeatMapColor(RNScalar value)
{
  // Compute color
  GLdouble r, g, b;
  if (value < 0.2) {
    value *= 4;
    r = 1;
    g = value;
    b = 0;
  }
  else if (value < 0.4) {
    value = (value - 0.25) * 4;
    r = 1 - value;
    g = 1;
    b = 0;
  }
  else if (value < 0.6) {
    value = (value - 0.5) * 4;
    r = 0;
    g = 1;
    b = value;
  }
  else if (value < 0.8) {
    value = (value - 0.75) * 4;
    r = 0;
    g = 1 - value;
    b = 1;
  }
  else {
    value = (value - 0.75) * 4;
    r = 0;
    g = 0;
    b = 1;
  }

  // Load color
  RNLoadRgb(r, g, b);
}



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static R3Vector
ComputeNormal(GSVScanline *scanline, int point_index)
{
  // For now
  if (load_points_from_images) return R3zero_vector;

  // Get convenient variables
  GSVScan *scan = scanline->Scan();
  int scanline_index = scanline->ScanIndex();
  R3Point p = scanline->PointPosition(point_index);

  // Find point on previous row
  R3Point pA = (point_index > 0) ? scanline->PointPosition(point_index-1) : p;

  // Find point on next row
  R3Point pB = (point_index < scanline->NPoints()-1) ? scanline->PointPosition(point_index+1) : p;

  // Find point on previous scanline
  R3Point p0 = p;
  GSVScanline *scanline0 = NULL;
  if (scanline_index > 0) {
    RNScalar closest_dd = FLT_MAX;
    scanline0 = scan->Scanline(scanline_index-1);
    for (int i = 0; i < scanline0->NPoints(); i++) {
      const R3Point& q = scanline0->PointPosition(i);
      RNLength dd = R3SquaredDistance(q, p);
      if (dd < closest_dd) { p0 = q; closest_dd = dd; }
    }
  }

  // Find point on next scanline
  R3Point p1 = p;
  GSVScanline *scanline1 = NULL;
  if (scanline_index < scan->NScanlines()-1) {
    RNScalar closest_dd = FLT_MAX;
    scanline1 = scan->Scanline(scanline_index-1);
    for (int i = 0; i < scanline1->NPoints(); i++) {
      const R3Point& q = scanline1->PointPosition(i);
      RNLength dd = R3SquaredDistance(q, p);
      if (dd < closest_dd) { p0 = q; closest_dd = dd; }
    }
  }

  // Compute normal with cross product
  R3Vector v01 = p1 - p0;
  v01.Normalize();
  R3Vector vAB = pB - pA;
  vAB.Normalize();
  R3Vector n = v01 % vAB;
  n.Normalize();

  // Flip to face laser
  R3Vector v = p - scanline->Pose().Viewpoint();
  if (n.Dot(v) > 0) n.Flip();

  // Return normal
  return n;
}



////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

static void
DrawFeature(GSVFeature *feature, RNBoolean transform)
{
  // Set transformation
  if (transform) {
    if (feature->image) PushTransformation(feature->image);
    else if (feature->scanline) PushTransformation(feature->scanline);
  }

  // Draw feature
  feature->Draw();

  // Reset transformation
  if (transform) PopTransformation();
}



static void 
DrawImage(GSVImage *image, int side, RNBoolean transform)
{
  // Persistent variables
  static GLuint texture_id[2] = { 0, 0 };
  static GSVImage *last_image[2] = { NULL, NULL };

  // Load texture
  if (image != last_image[side]) {
    last_image[side] = image;
    if (texture_id[side] > 0) glDeleteTextures(1, &texture_id[side]);
    texture_id[side] = 0;

    // Create texture
    R2Image *texture = ReadUndistortedImage(image);
    if (texture) {
      // Create texture
      glGenTextures(1, &texture_id[side]);
      glBindTexture(GL_TEXTURE_2D, texture_id[side]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      assert(texture_id[side] != 0);
      delete texture;
    }
  }

  // Draw textured polygon
  if (texture_id[side] > 0) {
    // Set transformation
    if (transform) PushTransformation(image);

    // Disable depth buffer 
    // glDisable(GL_DEPTH_TEST);
    // glDepthMask(FALSE);

    // Enable texture
    glBindTexture(GL_TEXTURE_2D, texture_id[side]);
    glEnable(GL_TEXTURE_2D);

    // Get camera parameters
    GSVCamera *camera = image->Camera();
    RNAngle xfov = 0.5 * camera->XFov();
    RNAngle yfov = 0.5 * camera->YFov();
    GSVPose pose = image->Pose();
    R3Point viewpoint = pose.Viewpoint();
    R3Vector towards = pose.Towards();
    R3Vector up = pose.Up();
    R3Vector right = pose.Right();

    // Draw textured polygon
    R3Point origin = viewpoint + towards * image_plane_distance;
    R3Vector dx = right * image_plane_distance * tan(xfov);
    R3Vector dy = up * image_plane_distance * tan(yfov);
    R3Point ur = origin + dx + dy;
    R3Point lr = origin + dx - dy;
    R3Point ul = origin - dx + dy;
    R3Point ll = origin - dx - dy;
    glColor3d(1,1,1);
    glBegin(GL_POLYGON);
    glTexCoord2d(0,0);
    R3LoadPoint(ll);
    glTexCoord2d(1,0);
    R3LoadPoint(lr);
    glTexCoord2d(1, 1);
    R3LoadPoint(ur);
    glTexCoord2d(0,1);
    R3LoadPoint(ul);
    glEnd();

    // Disable texture
    glDisable(GL_TEXTURE_2D);

    // Enable depth buffer
    // glEnable(GL_DEPTH_TEST);
    // glDepthMask(TRUE);

    // Reset transformation
    if (transform) PopTransformation();
  }
}



static void 
DrawFeatures(int side, RNBoolean transform, int color_scheme)
{
  // Draw features
  for (int i = 0; i < optimization->NFeatures(); i++) {
    GSVFeature *feature = optimization->Feature(i);

    // Check if should draw feature
    if (feature->image && selected_image[side] && (feature->image != selected_image[side])) continue;

    // Draw feature
    LoadFeatureColor(feature, side, color_scheme, i);
    DrawFeature(feature, transform);
  }
}



static void 
DrawCorrespondences(int side, RNBoolean transform, int color_scheme)
{
  // Draw correspondences
  for (int i = 0; i < optimization->NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);

    // Check if should draw correspondence
    if (selected_image[side]) {
      if (feature0->image && (feature0->image != selected_image[side]) && 
        feature1->image && (feature1->image != selected_image[side])) continue;
    }
  if (selected_image[1-side]) {
    if (feature0->image && (feature0->image != selected_image[1-side]) && 
      feature1->image && (feature1->image != selected_image[1-side])) continue;
  }

    // Draw features
LoadFeatureColor(feature0, side, color_scheme, 0);
DrawFeature(feature0, transform);
LoadFeatureColor(feature1, side, color_scheme, 0);
DrawFeature(feature1, transform);

    // Draw line between features
glBegin(GL_LINES);
LoadFeatureColor(feature0, side, color_scheme, 0);
R3LoadPoint(optimization->FeaturePosition(feature0, transform));
LoadFeatureColor(feature1, side, color_scheme, 0);
R3LoadPoint(optimization->FeaturePosition(feature1, transform));
glEnd();
}
}


static void 
DrawScans(int side, RNBoolean transform, int color_scheme)
{
  // Draw scans
  glBegin(GL_POINTS);
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      if (show_all_runs || ir == run_number) {
       GSVRun *run = scene->Run(ir);
       for (int is = 0; is < run->NSegments(); is++) {
         if (show_all_segments || is == segment_number) {
           GSVSegment *segment = run->Segment(is);
           LoadSegmentColor(segment, color_scheme);
           for (int ia = 0; ia < segment->NScans(); ia++) {
             GSVScan *scan = segment->Scan(ia);
             if (scan->NScanlines() == 0) continue;
             ScanData *scan_data = (ScanData *) scan->Data();
             if (!scan_data) continue;
             for (int ie = 0; ie < scan->NScanlines(); ie++) {
              GSVScanline *scanline = scan->Scanline(ie);
              double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
              if (timestamp_delta > timestamp_radius) continue;
              const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
              for (int ik = 0; ik < scanline->NPoints(); ik++) {
                R3Point position = scanline->PointPosition(ik);
                if (transform) position.Transform(scanline_transformation);
                R3LoadPoint(position);
              }
            }
          }
        }
      }
    }
  }
}
else if (color_scheme == VIEWPOINT_DEPTH_COLOR_SCHEME) {
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
     GSVRun *run = scene->Run(ir);
     for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           if (scan->NScanlines() == 0) continue;
           ScanData *scan_data = (ScanData *) scan->Data();
           if (!scan_data) continue;
           for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            R3Point viewpoint = scanline->Pose().Viewpoint();
            R3Vector towards = scanline->Pose().Towards();
            if (transform) viewpoint.Transform(scanline_transformation);
            if (transform) towards.Transform(scanline_transformation);
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              if (transform) position.Transform(scanline_transformation);
              R3Vector vector = position - viewpoint;
              RNScalar depth = vector.Dot(towards);
              RNScalar value = depth / image_plane_distance;
              LoadHeatMapColor(value);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
}
}
else if (color_scheme == HEIGHT_COLOR_SCHEME) {
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
     GSVRun *run = scene->Run(ir);
     for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           if (scan->NScanlines() == 0) continue;
           ScanData *scan_data = (ScanData *) scan->Data();
           if (!scan_data) continue;
           for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            RNCoord ground_z = scanline->EstimatedGroundZ();
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              RNScalar height = position.Z() - ground_z;
              if (transform) position.Transform(scanline_transformation);
              RNScalar value = 2 * height / image_plane_distance;
              LoadHeatMapColor(value);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
}
}
else if (color_scheme == PICK_COLOR_SCHEME) {
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
     GSVRun *run = scene->Run(ir);
     for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           if (scan->NScanlines() == 0) continue;
           ScanData *scan_data = (ScanData *) scan->Data();
           if (!scan_data) continue;
           for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
            if (!scanline_data) continue;
            int scanline_index = scanline_data->index;
            unsigned char rgba[4];
            rgba[0] = (scanline_index >> 16) & 0xFF;
            rgba[1] = (scanline_index >> 8) & 0xFF;
            rgba[2] = scanline_index & 0xFF;
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              if (transform) position.Transform(scanline_transformation);
              rgba[3] = (ik + 1) & 0xFF;
              glColor4ubv(rgba);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
}
}
glEnd();
}



static void 
DrawMeshes(int side, RNBoolean transform, int color_scheme)
{
  // Colors
  const int max_mesh_colors = 12;
  const RNRgb mesh_colors[max_mesh_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),
    
  };

  // Set OpenGL modes
  glEnable(GL_LIGHTING);
  static GLfloat material[] = { 0.8, 0.8, 0.8, 1.0 };

  // Draw meshes
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         RNRgb color = mesh_colors[(7*ir + is) % max_mesh_colors];
         material[0] = color[0]; material[1] = color[1]; material[2] = color[2];
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material); 
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           ScanData *scan_data = (ScanData *) scan->Data();
           if (!scan_data) continue;
           if (scan->NScanlines() == 0) continue;
           if (!scan_data->mesh) continue;
           scan_data->mesh->Draw();
         }
       }
     }
   }
 }

  // Reset OpenGL modes
 glDisable(GL_LIGHTING);
}



static void 
DrawLasers(int side, RNBoolean transform, int color_scheme)
{
    // Scanline drawing parameters
  static int scanline_step = 0;
  if (scanline_step == 0) {
    const int max_scanlines = 32 * 1024;
    scanline_step = 1 + scene->NScanlines() / max_scanlines;
  }
  
    // Draw pose for every laser
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           for (int ie = 0; ie < scan->NScanlines(); ie += scanline_step) {
             GSVScanline *scanline = scan->Scanline(ie);
             GSVPose pose = scanline->Pose();
             R3Point position = pose.Viewpoint();
             R3Vector towards = pose.Towards();
             R3Vector up = pose.Up();
             if (transform) {
              R3Affine transformation = optimization->OptimizedTransformation(scanline);
              position.Transform(transformation);
              towards.Transform(transformation);
              up.Transform(transformation);
            }
            if (scanline == selected_scanline[side]) glColor3d(1.0, 1.0, 1.0);
            else LoadScanlineColor(scanline, side, color_scheme);
            R3Span(position, position + towards).Draw();
            R3Span(position, position + 0.5*up).Draw();
          }
        }
      }
    }
  }
}
}



static void 
DrawCameras(int side, RNBoolean transform, int color_scheme)
{
  // Draw camera pose for every image
  int image_index = 0;
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ip = 0; ip < segment->NPanoramas(); ip++) {
           GSVPanorama *panorama = segment->Panorama(ip);
           for (int ii = 0; ii < panorama->NImages(); ii++) {
             GSVImage *image = panorama->Image(ii);
             GSVPose pose = image->Pose();
             R3Point position = pose.Viewpoint();
             R3Vector towards = pose.Towards();
             R3Vector up = pose.Up();
             if (transform) {
              R3Affine transformation = optimization->OptimizedTransformation(image);
              position.Transform(transformation);
              towards.Transform(transformation);
              up.Transform(transformation);
            }
            LoadImageColor(image, side, color_scheme);
            R3Span(position, position + towards).Draw();
            R3Span(position, position + 0.5 * up).Draw();
            image_index++;
          }
        }
      }
    }
  }
}
}



static void 
DrawPaths(int side, RNBoolean transform, int color_scheme)
{
  // Draw path for every segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {

      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         SegmentData *segment_data = (SegmentData *) segment->Data();
         if (!segment_data) continue;
         GSVPath *path = segment_data->path;
         LoadSegmentColor(segment, color_scheme);

	  // Draw vertices
         glPointSize(3);
         glBegin(GL_POINTS);
         for (int i = 0; i < path->NVertices(); i++) {
           R3Point point = path->VertexPosition(i);
           if (transform) point += path->VertexTranslation(i);
           R3LoadPoint(point);
         }
         glEnd();
         glPointSize(1);

	  // Draw spline curve
         glBegin(GL_LINE_STRIP);
         for (RNScalar u = 0; u <= path->VertexParameter(path->NVertices()-1); u += 0.1) {
           R3Point point = path->spline->PointPosition(u);
           if (transform) point += path->Translation(u);
           R3LoadPoint(point);
         }
         glEnd();
       }
     }
   }
 }
}



static void 
DrawSegments(int side, RNBoolean transform, int color_scheme)
{
  // Colors
  const int max_segment_colors = 12;
  const RNRgb segment_colors[max_segment_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),

  };

  // Draw segments
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
       if (show_all_segments || is == segment_number) {
         GSVSegment *segment = run->Segment(is);
         for (int ia = 0; ia < segment->NScans(); ia++) {
           GSVScan *scan = segment->Scan(ia);
           ScanData *scan_data = (ScanData *) scan->Data();
           if (!scan_data) continue;
           if (scan->NScanlines() == 0) continue;
           for (int i = 0; i < scan_data->da_position_x_grid.XResolution(); i++) {
             for (int j = 0; j < scan_data->da_position_x_grid.YResolution(); j++) {
              RNScalar x = scan_data->da_position_x_grid.GridValue(i, j);
              if (x == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar y = scan_data->da_position_y_grid.GridValue(i, j);
              if (y == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar z = scan_data->da_position_z_grid.GridValue(i, j);
              if (z == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar scanline_index_value = scan_data->da_scanline_grid.GridValue(i, j);
              if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar segment_index_value = scan_data->da_plane_segment_id_grid.GridValue(i, j);
              if (segment_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              int scanline_index = (int) (scanline_index_value + 0.5);
              int segment_index = (int) (segment_index_value + 0.5);
              GSVScanline *scanline = scan->Scanline(scanline_index);
              const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
              R3Point position(x, y, z);
              if (transform) position.Transform(scanline_transformation);
              RNLoadRgb(segment_colors[segment_index % max_segment_colors]);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
}
glEnd();
}



static void 
DrawImages(int side, RNBoolean transform, int color_scheme)
{
  // Check selected image
  if (!selected_image[side]) return;

  // Draw selected image
  DrawImage(selected_image[side], side, transform);
}



static void 
DrawBBoxes(int side, RNBoolean transform, int color_scheme)
{
  // Draw bbox for every segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    if (show_all_runs || ir == run_number) {
      if (ir != run_count) continue;
      GSVRun *run = scene->Run(ir);
      //run->BBox().Outline();
      for (int is = 0; is < run->NSegments(); is++) {
       	if (show_all_segments || is == segment_number) {
       	  GSVSegment *segment = run->Segment(is);
       	  //LoadSegmentColor(segment, color_scheme);
          if (is != segment_count) continue;
          printf("%s %d\n", run->Name(), is);

       	  segment->BBox().Outline();
       	}
       }
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Pick functions
////////////////////////////////////////////////////////////////////////

static int 
Pick(GSVPoseOptimization *optimization, int x, int y, 
  RNBoolean pick_features = TRUE, RNBoolean pick_cameras = TRUE, RNBoolean pick_lasers = TRUE, RNBoolean pick_scans = TRUE,
  GSVFeature **picked_feature = NULL, GSVImage **picked_image = NULL, GSVScanline **picked_scanline = NULL, int *picked_point_index = NULL, 
  R3Point *picked_position = NULL)
{
  // How close the cursor has to be to a point (in pixels)
  int pick_tolerance = 10;

  // Clear window 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set viewing transformation
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);
  int viewport_width = GLUTwindow_width/(split_screen+1);
  glViewport(side*viewport_width, 0, viewport_width, GLUTwindow_height);
  viewer[side]->Camera().Load();

  // Set OpenGL stuff
  glLineWidth(2 * pick_tolerance);
  glPointSize(2 * pick_tolerance);    

  // Draw everything
  if (pick_features && (show_features[side] || show_correspondences[side])) DrawFeatures(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_cameras && show_cameras[side]) DrawCameras(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_lasers && show_lasers[side]) DrawLasers(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_scans && show_scans[side]) DrawScans(side, show_optimization[side], PICK_COLOR_SCHEME);

  // Reset OpenGL stuff
  glPointSize(1);
  glLineWidth(1);
  glFinish();

  // Read color buffer at cursor position
  unsigned char rgba[4];
  glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
  if (rgba[3] == 0) return 0;

  // Return scanline
  int r = rgba[0] & 0xFF;
  int g = rgba[1] & 0xFF;
  int b = rgba[2] & 0xFF;
  int a = rgba[3] & 0xFF;

  // Determine pick type
  int pick_type = 0;
  GSVScanline *scanline = NULL;
  if (a == 0xFB) {
    // Picked feature
    int feature_index = (r << 16) | (g << 8) | b;
    if (feature_index < 0) return 0;
    if (feature_index >= optimization->NFeatures()) return 0;
    GSVFeature *feature = optimization->Feature(feature_index);
    if (picked_feature) *picked_feature = feature;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = feature->scanline;
    if (picked_point_index) *picked_point_index = feature->scan_point_index;
    scanline = feature->scanline;
    pick_type = 1;
  }
  else if (a == 0xFC) {
    // Picked camera
    int image_index = (r << 16) | (g << 8) | b;
    if (image_index < 0) return 0;
    if (image_index > optimization->images.NEntries()) return 0;
    GSVImage *image = optimization->images.Kth(image_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = image;
    if (picked_scanline) *picked_scanline = NULL;
    if (picked_point_index) *picked_point_index = -1;
    pick_type = 2;
  }
  else if (a == 0xFD) {
    // Picked laser
    int scanline_index = (r << 16) | (g << 8) | b;
    if (scanline_index < 0) return 0;
    if (scanline_index > optimization->scanlines.NEntries()) return 0;
    scanline = optimization->scanlines.Kth(scanline_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = scanline;
    if (picked_point_index) *picked_point_index = -1;
    pick_type = 3;
  }
  else {
    // Picked scan point
    pick_type = 1;
    int scanline_index = (r << 16) | (g << 8) | b;
    if (scanline_index < 0) return 0;
    if (scanline_index > optimization->scanlines.NEntries()) return 0;
    scanline = optimization->scanlines.Kth(scanline_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = scanline;
    if (picked_point_index) *picked_point_index = a - 1;
    pick_type = 4;
  }

  // Return position
  if (picked_position) {
    // Find hit position
    GLfloat depth;
    GLdouble p[3];
    GLint viewport[4];
    GLdouble modelview_matrix[16];
    GLdouble projection_matrix[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    gluUnProject(x, y, depth, modelview_matrix, projection_matrix, viewport, &(p[0]), &(p[1]), &(p[2]));
    R3Point position(p[0], p[1], p[2]);
    *picked_position = position;
  }
  
  // Return pick type
  return pick_type;
}



////////////////////////////////////////////////////////////////////////
// GLUT callback functions
////////////////////////////////////////////////////////////////////////

void AtExit(void)
{
  // Write correspondences
  if (output_correspondences_name) {
    WriteCorrespondences(optimization, output_correspondences_name);
  }

  // Write transformation
  if (output_transformations_name) {
    if (!WriteTransformations(optimization, output_transformations_name)) exit(-1);
  }

  // Write error plot
  if (output_error_plot_name) {
    WriteErrorPlot(optimization, output_error_plot_name);
  }

  // Write scene
  if (output_scene_name) {
    WriteScene(scene, optimization, output_scene_name);
  }
}



void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void GLUTRedraw(void)
{
  // Clear window 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw both sides of screen
  for (int side = 0; side <= split_screen; side++) {
    // Update viewer to match selected camera
    if (viewpoint_scheme[side] && selected_image[side]) {
      GSVCamera *camera = selected_image[side]->Camera();
      if (camera) {
        RNAngle xfov = 0.5 * camera->XFov();
        RNAngle yfov = 0.5 * camera->YFov();
        GSVPose pose = selected_image[side]->Pose();
        R3Point viewpoint = pose.Viewpoint();
        R3Vector towards = pose.Towards();
        R3Vector up = pose.Up();
        if (show_optimization[side]) {
          R3Affine transformation = optimization->OptimizedTransformation(selected_image[side]);
          viewpoint.Transform(transformation);
          towards.Transform(transformation);
          up.Transform(transformation);
        }
        center_point[side] = viewpoint + towards;
        R3Camera c(viewpoint, towards, up, xfov, yfov, 0.001, 1000000);
        viewer[side]->SetCamera(c);
      }
    }

    // Set viewport
    if (!split_screen) glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);    
    else glViewport(side*GLUTwindow_width/2, 0, GLUTwindow_width/2, GLUTwindow_height);
    
    // Set viewing transformation
    viewer[side]->Camera().Load();
    glPointSize(point_size[side]);

    // Draw correspondences
    if (show_correspondences[side]) {
      DrawCorrespondences(side, show_optimization[side], color_scheme[side]);
    }

    // Draw features
    if (show_features[side]) {
      // Draw features
      DrawFeatures(side, show_optimization[side], color_scheme[side]);
    }

    // Draw selected feature
    if (selected_feature[side]) {
      RNLoadRgb(1.0, 1.0, 1.0);
      DrawFeature(selected_feature[side], show_optimization[side]);
    }
    
    // Draw scans
    if (show_scans[side]) {
      DrawScans(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw meshes
    if (show_meshes[side]) {
      DrawMeshes(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw paths
    if (show_paths[side]) {
      DrawPaths(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw cameras
    if (show_cameras[side]) {
      DrawCameras(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw lasers
    if (show_lasers[side]) {
      DrawLasers(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw images
    if (show_images[side]) {
      DrawImages(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw segments
    if (show_segments[side]) {
      DrawSegments(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw bounding box
    if (show_bboxes[side]) {
      DrawBBoxes(side, show_optimization[side], color_scheme[side]);
    }
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  if (split_screen) {
    viewer[0]->ResizeViewport(0, 0, w/2, h);
    viewer[1]->ResizeViewport(w/2, 0, w/2, h);
  }
  else {
    viewer[0]->ResizeViewport(0, 0, w, h);
    viewer[1]->ResizeViewport(0, 0, w, h);
  }

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}



void GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute mouse movement
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];
  
  // World in hand navigation 
  if (GLUTbutton[0]) viewer[GLUTside]->RotateWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  else if (GLUTbutton[1]) viewer[GLUTside]->ScaleWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  else if (GLUTbutton[2]) viewer[GLUTside]->TranslateWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Check for drag
  static int last_mouse_down[2] = { -999, -999 };
  RNBoolean drag = (fabs(x - last_mouse_down[0]) > 10) | (fabs(y - last_mouse_down[1]) > 10);
  last_mouse_down[0] = (state == GLUT_DOWN) ? x : -999;
  last_mouse_down[1] = (state == GLUT_DOWN) ? y : -999;

  // Process mouse button event
  if (!drag && (state == GLUT_UP)) {
    if (button == GLUT_LEFT_BUTTON) {
      // Check for double click
      static RNBoolean double_click = FALSE;
      static RNTime last_mouse_down_time;
      double_click = !double_click && (last_mouse_down_time.Elapsed() < 1);
      last_mouse_down_time.Read();

      if (double_click) {
        // Toggle display of selected segment
        R3Point position;
        GSVScanline *scanline = NULL;
        if (Pick(optimization, x, y, FALSE, FALSE, show_lasers[side], show_scans[side], NULL, NULL, &scanline, NULL, &position)) {
          center_point [side] = position;
          if (scanline) {
            RNScalar timestamp = scanline->Timestamp();
            GSVScan *scan = scanline->Scan();
            GSVSegment *segment = scan->Segment();
            ToggleScanData(segment, timestamp);
            glutPostRedisplay();
          }
        }
      }
      else {
        // Set center point
       bool split = false;
       if (viewpoint_scheme[side] && selected_image[side]) {
	  // Add image feature
         int sx = (side == 0) ? x : x - (GLUTwindow_width/2);
         double ix = (double) (sx * selected_image[side]->Width()) / (GLUTwindow_width/2);
         double iy = (double) (y * selected_image[side]->Height()) / GLUTwindow_height;
         GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, 
          selected_image[side], R2Point(ix, iy), R2zero_vector, 1.0, 1.0);
         optimization->InsertFeature(feature);
         last_feature = feature;
         last_correspondence = NULL;
         split = true;
         if (split_screen && selected_feature[1-side]) {
           GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1000);
           optimization->InsertCorrespondence(correspondence);
           last_correspondence = correspondence;
         }
         else if (selected_feature[side]) {
           GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1000);
           optimization->InsertCorrespondence(correspondence);
           last_correspondence = correspondence;
         }
       }
       if (!split) {
         R3Point position;
         if (Pick(optimization, x, y, show_features[side], show_cameras[side], show_lasers[side], show_scans[side], NULL, NULL, NULL, NULL, &position)) {
           center_point[side] = position;
           glutPostRedisplay();
           printf("%f %f %f\n", position.X(), position.Y(), position.Z());
         }
       }
       else {
         glutPostRedisplay();
       }
     }
   }
 }

  // Remember button state 
 int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
 GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
 GLUTmodifiers = glutGetModifiers();

  // Remember mouse position 
 GLUTmouse[0] = x;
 GLUTmouse[1] = y;

  // Remember side 
 GLUTside = side;
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Process keyboard button event 
  switch (key) {
    case 'A':
    case 'a': 
    case 'N':
    case 'n': 
    if (1) {
      // Add scan feature
      R3Point position;
      int point_index = -1;
      GSVScanline *scanline = NULL;
      if (Pick(optimization, x, y, FALSE, FALSE, FALSE, show_scans[side], NULL, NULL, &scanline, &point_index, &position)) {
        if (point_index >= 0) {
          center_point[side] = position;
          R3Vector normal = ComputeNormal(scanline, point_index);
          if (show_optimization[side]) position.InverseTransform(optimization->OptimizedTransformation(scanline));
          int feature_type = ((key == 'a') || (key == 'A')) ? GSV_SCAN_POINT_FEATURE_TYPE : GSV_SCAN_PLANE_FEATURE_TYPE;
          GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, R3zero_vector, normal, 1.0, 1.0);
          optimization->InsertFeature(feature);
          last_feature = feature;
          last_correspondence = NULL;
          if ((key == 'A') || (key == 'N')) {
            if (split_screen && selected_feature[1-side]) {
              GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1000);
              optimization->InsertCorrespondence(correspondence);
              last_correspondence = correspondence;
            }
            else if (selected_feature[side]) {
              GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1000);
              optimization->InsertCorrespondence(correspondence);
              last_correspondence = correspondence;
            }
          }
          selected_feature[side] = feature;
        }
      }
    }
    break; 

    case 'X':
    case 'x': 
    if (selected_image[side]) {
      // Add image feature
      int sx = (side == 0) ? x : x - (GLUTwindow_width/2);
      double ix = (double) (sx * selected_image[side]->Width()) / (GLUTwindow_width/2);
      double iy = (double) (y * selected_image[side]->Height()) / GLUTwindow_height;
      GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, 
        selected_image[side], R2Point(ix, iy), R2zero_vector, 1.0, 1.0);
      optimization->InsertFeature(feature);
      last_feature = feature;
      last_correspondence = NULL;
      if (key == 'X') {
        if (split_screen && selected_feature[1-side]) {
          GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1000);
          optimization->InsertCorrespondence(correspondence);
          last_correspondence = correspondence;
        }
        else if (selected_feature[side]) {
          GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1000);
          optimization->InsertCorrespondence(correspondence);
          last_correspondence = correspondence;
        }
      }
      selected_feature[side] = feature;
    }
    break; 


    case 'R':
    case 'r':
    run_count++;
    break;

    case 'T':
    case 't':
    run_count--;
    break;

    case 'Y':
    case 'y':
    segment_count++;
    break;

    case 'U':
    case 'u':
    segment_count--;
    break;

    case 'B':
    case 'b':
    show_bboxes[side] = !show_bboxes[side];
    break;

    case 'C':
    case 'c':
    show_cameras[side] = !show_cameras[side];
    break;

    case 'D':
    case 'd':
    if (color_scheme[side] == SEGMENT_COLOR_SCHEME) color_scheme[side] = VIEWPOINT_DEPTH_COLOR_SCHEME;
    else if (color_scheme[side] == VIEWPOINT_DEPTH_COLOR_SCHEME) color_scheme[side] = HEIGHT_COLOR_SCHEME;
    else color_scheme[side] = SEGMENT_COLOR_SCHEME;
    break;

    case 'E':
    case 'e':
    WriteErrorPlot(optimization);
    break;

    case 'F':
    case 'f':
    show_features[side] = !show_features[side];
    break;

  // print feature id
    case 'G':
    case 'g':
    if (show_features[side] || show_correspondences[side]) {
      GSVFeature *feature = NULL;
      if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature)) {
        if (feature) {
         printf("Selected Feature: %d\n", feature->index);
       }
     }
   }
   break;

  // output intersection coordinate
   case 'H':
   case 'h':
   if (1) {
    FILE *fp = fopen("intersections.txt", "a");

      // Draw both sides of screen
    for (int i = 0; i <= split_screen; i++) {
	// Update viewer to match selected camera
     if (viewpoint_scheme[i] && selected_image[i]) {
       GSVCamera *camera = selected_image[i]->Camera();
       if (camera) {
         GSVPose pose = selected_image[i]->Pose();
         R3Point viewpoint = pose.Viewpoint();
         fprintf(fp, "Intersection %d - %f %f %f\n", total_intersections, viewpoint.X(), viewpoint.Y(), viewpoint.Z());
         total_intersections++;
       }
     }
   }
   fclose(fp);
   break;
 }

 case 'I':
 case 'i':
 show_images[side] = !show_images[side];
 break;

 case 'L':
 case 'l':
 show_lasers[side] = !show_lasers[side];
 break;

 case 'M':
 case 'm':
 show_correspondences[side] = !show_correspondences[side];
 break;

 case 'O':
 case 'o':
 show_optimization[side] = !show_optimization[side];
 break;

 case 'P':
 case 'p':
 show_scans[side] = !show_scans[side];
 break;

 case 'S':
 case 's':
 show_paths[side] = !show_paths[side];
 break;

 case 'V':
 case 'v':
 viewpoint_scheme[side] = !viewpoint_scheme[side];
 break;

 case 'W':
 case 'w':
 split_screen = 1 - split_screen;
 GLUTResize(GLUTwindow_width, GLUTwindow_height);
 break;

 case '%':
 solver = (solver + 1) % RN_POLYNOMIAL_NUM_SOLVERS;
 printf("Solver = %d\n", solver);
 break;

 case '^':
 show_segments[side] = !show_segments[side];
 break;

 case '&':
 show_meshes[side] = !show_meshes[side];
 break;

 case '#':
    // optimization->CreateImageToScanCorrespondences();
 break;

 case '+':
 point_size[side] += 1;
 break;

 case '_':
 point_size[side] -= 1;
 break;

 case '=':
 if (image_plane_distance < 1000) image_plane_distance *= 1.1;
 break;

 case '-':
 if (image_plane_distance > 0.01) image_plane_distance *= 0.9;
 break;

  case 1:  // ctrl-A
    // Update poses
  Solve(optimization, TRUE, TRUE,  TRUE, TRUE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 2:  // ctrl-B
    // Update all translations (including camera and laser offsets)
  Solve(optimization, TRUE, FALSE,  TRUE, TRUE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 3:  // ctrl-C
    // Update camera translations only
  Solve(optimization, TRUE, FALSE,  FALSE, TRUE, FALSE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 4:  // ctrl-D
    // Update pixel depths only
  Solve(optimization, FALSE, FALSE,  FALSE, FALSE, FALSE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 5: // ctrl-E
    // Clear pose transformations
  optimization->ClearOptimizationVariables();
  optimization->UpdatePixelPositionsFromCorrespondenceSolutions();
  break;

  case 6:  // ctrl-F
    // Update pose translations and rotations
  Solve(optimization, TRUE, TRUE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 12:  // ctrl-L
    // Update laser translations only
  Solve(optimization, TRUE, FALSE,  TRUE, FALSE, FALSE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 15:  // ctrl-O
    // Randomize pose translations
  RandomizePoseTranslations(optimization);
  break;

  case 18:  // ctrl-R
    // Update pose rotations
  Solve(optimization, FALSE, TRUE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 19:  // ctrl-S
  if (output_correspondences_name) {
    WriteCorrespondences(optimization, output_correspondences_name);
  }
  break;

  case 20:  // ctrl-T
    // Update pose translations
  Solve(optimization, TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity, solver);
  break;

  case 26:  // ctrl-Z
  if (last_correspondence) {
      // Remove last correspondence
    printf("HEREA\n");
    optimization->RemoveCorrespondence(last_correspondence);
    delete last_correspondence;
    last_correspondence = NULL;
  }
  if (last_feature) {
      // Remove last feature
    printf("HEREB\n");
    if (last_feature == selected_feature[0]) selected_feature[0] = NULL;
    if (last_feature == selected_feature[1]) selected_feature[1] = NULL;
    optimization->RemoveFeature(last_feature);
    delete last_feature;
    last_feature = NULL;
  }
  break;

  case 17: // ctrl-Q
    // Quit
  GLUTStop();
  break;

  case ' ':
  if (1) {
      // Select segment and unselect all others
    R3Point position;
    GSVImage *image = NULL;
    GSVScanline *scanline = NULL;
    if (Pick(optimization, x, y, FALSE, show_cameras[side], show_lasers[side], show_scans[side], NULL, &image, &scanline, NULL, &position)) {
      center_point[side] = position;
      GSVSegment *segment = NULL;
      RNScalar timestamp = -1;
      if (scanline) {
        selected_scanline[side] = scanline;
        segment = scanline->Scan()->Segment();
        timestamp = scanline->Timestamp();
      }
      else if (image) {
        selected_image[side] = image;
        segment = image->Tapestry()->Segment();
        timestamp = image->Timestamp();
      }
      if (segment) {
        if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) ToggleScanData(segment, timestamp);
        else SelectScanData(segment, timestamp);
      }
    }
  }
  break;

  case 27: // ESC
    // Select no segments
  selected_feature[0] = NULL;
  selected_feature[1] = NULL;
  selected_scanline[0] = NULL;
  selected_scanline[1] = NULL;
  selected_image[0] = NULL;
  selected_image[1] = NULL;
  SelectScanData(NULL, -1);
  break;
}

  // Remember mouse position 
GLUTmouse[0] = x;
GLUTmouse[1] = y;

  // Remember modifiers 
GLUTmodifiers = glutGetModifiers();

  // Remember side 
GLUTside = side;

  // Redraw
glutPostRedisplay();  
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Process keyboard button event 
  switch (key) {
    case GLUT_KEY_F1:
    // Select feature
    if (show_features[side] || show_correspondences[side]) {
      selected_feature[side] = NULL;
      if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &selected_feature[side])) {
        if (selected_feature[side]) {
          if ((selected_feature[side]->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) ||
            (selected_feature[side]->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) ||
            (selected_feature[side]->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
            GSVScanline *scanline = selected_feature[side]->scanline;
          GSVScan *scan = (scanline) ? scanline->Scan() : NULL;
          GSVSegment *segment = (scan) ? scan->Segment() : NULL;
          GSVRun *run = (segment) ? segment->Run() : NULL;
          int scanline_index = (scanline) ? scanline->ScanIndex() : -1;
          int scan_index = (scan) ? scan->SegmentIndex() : -1;
          int segment_index = (segment) ? segment->RunIndex() : -1;
          const char *run_name = (run) ? run->Name() : "None";
          R3Point scan_position = optimization->FeaturePosition(selected_feature[side], FALSE); 
          R3Vector scan_direction = optimization->FeatureNormal(selected_feature[side], FALSE); 
          GSVDescriptor *descriptor = &selected_feature[side]->descriptor;
          printf("%6d   %30s %2d %2d %6d %3d   %8.6f   %12.6f %12.6f %12.6f  %8.6f %8.6f %8.6f\n", 
           selected_feature[side]->index, 
           run_name, segment_index, scan_index, scanline_index, 
           selected_feature[side]->scan_point_index, selected_feature[side]->score,
           scan_position.X(), scan_position.Y(), scan_position.Z(), 
           scan_direction.X(), scan_direction.Y(), scan_direction.Z());
          for (int i = 0; i < descriptor->nvalues; i++) 
            printf("%6.3f ", descriptor->values[i]);
          printf("\n");
        }
        else if ((selected_feature[side]->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
         (selected_feature[side]->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) {
          GSVImage *image = selected_feature[side]->image;
        GSVPanorama *panorama = (image) ? image->Panorama() : NULL;
        GSVSegment *segment = (panorama) ? panorama->Segment() : NULL;
        GSVRun *run = (segment) ? segment->Run() : NULL;
        int image_index = (image) ? image->PanoramaIndex() : -1;
        int panorama_index = (panorama) ? panorama->SegmentIndex() : -1;
        int segment_index = (segment) ? segment->RunIndex() : -1;
        const char *run_name = (run) ? run->Name() : "None";
        GSVDescriptor *descriptor = &selected_feature[side]->descriptor;
        printf("%6d   %30s %2d %6d %2d   %8.6f   %12.6f %12.6f   %8.6f %8.6f\n", 
         selected_feature[side]->index, 
         run_name, segment_index,  panorama_index, image_index, selected_feature[side]->score,
         selected_feature[side]->image_position.X(), selected_feature[side]->image_position.Y(), 
         selected_feature[side]->image_direction.X(), selected_feature[side]->image_direction.Y());
        for (int i = 0; i < descriptor->nvalues; i++) 
          printf("%6.3f ", descriptor->values[i]);
        printf("\n");
      }
    }
  }
}
break; 

case GLUT_KEY_F2:
    // Create correspondence to selected feature 
if (show_features[side] || show_correspondences[side]) {
  GSVFeature *feature = NULL;
  if (split_screen && selected_feature[1-side]) {
    if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature)) {
      if (feature != selected_feature[1-side]) {
        GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1000);
        optimization->InsertCorrespondence(correspondence);
        last_correspondence = correspondence;
        last_feature = NULL;
      }
    }
  }
  else if (selected_feature[side]) {
    if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature)) {
      if (feature != selected_feature[1-side]) {
        GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1000);
        optimization->InsertCorrespondence(correspondence);
        last_correspondence = correspondence;
        last_feature = NULL;
      }
    }
  }
}
break; 

case GLUT_KEY_F3:
    // Delete all correspondences with selected feature
if (show_features[side] || show_correspondences[side]) {
  GSVFeature *feature = NULL;
  if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature)) {
    if (feature) {
      RNArray<GSVFeatureCorrespondence *> correspondences;
      for (int i = 0; i < optimization->NCorrespondences(); i++)
        correspondences.Insert(optimization->Correspondence(i));
      for (int i = 0; i < correspondences.NEntries(); i++) {
        GSVFeatureCorrespondence *correspondence = correspondences.Kth(i);
        if ((correspondence->Feature(0) == feature) || (correspondence->Feature(1) == feature)) {
          optimization->RemoveCorrespondence(correspondence);
          delete correspondence;
        }
      }
    }
  }
}
break; 

case GLUT_KEY_F5:
rigidity = 0.75 * rigidity;
printf("Rigidity = %g\n", rigidity);
break;

case GLUT_KEY_F6:
rigidity = 1.5 * rigidity;
printf("Rigidity = %g\n", rigidity);
break;

case GLUT_KEY_F7:
optimization->max_intracorrespondence_euclidean_distance *= 0.75;
printf("Max distance = %g\n", optimization->max_intracorrespondence_euclidean_distance);
break;

case GLUT_KEY_F8:
optimization->max_intracorrespondence_euclidean_distance *= 1.5;
printf("Max distance = %g\n", optimization->max_intracorrespondence_euclidean_distance);
break;

case GLUT_KEY_F9: {
  static int correspondence_index = 0;
  if (correspondence_index < optimization->NCorrespondences()) {
    split_screen = 1;
    GSVFeatureCorrespondence *correspondence = optimization->Correspondence(correspondence_index);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    selected_feature[0] = feature0;
    selected_feature[1] = feature1;
    selected_image[0] = feature0->image;
    selected_image[1] = feature1->image;
    GLUTResize(GLUTwindow_width, GLUTwindow_height);
    GLUTKeyboard('v', GLUTwindow_width/4, GLUTwindow_width/2);
    GLUTKeyboard('v', 3*GLUTwindow_width/4, GLUTwindow_width/2);
    correspondence_index = (correspondence_index < optimization->NCorrespondences()-1) ? correspondence_index + 1 : 0;
  }
  break; }

  case GLUT_KEY_F10:
  optimization->correspondences.correspondences.Truncate(num_initial_correspondences);
  optimization->CreateMutuallyClosestCorrespondences();
  printf("%d\n", optimization->NCorrespondences());
  break;

  case GLUT_KEY_F11:
  optimization->correspondences.correspondences.Truncate(num_initial_correspondences);
  optimization->CreateMutuallyClosestCorrespondences();
  optimization->Solve(TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity);
  printf("%g %g %d %g\n", optimization->max_intracorrespondence_euclidean_distance, rigidity, 
    optimization->NCorrespondences(), optimization->Score());
  break;

  case GLUT_KEY_F12: {
    static int first = 1;
    if (first) { optimization->max_intracorrespondence_euclidean_distance = 20; rigidity = 1; first = 0; }
    optimization->correspondences.correspondences.Truncate(num_initial_correspondences);
    optimization->CreateMutuallyClosestCorrespondences();
    optimization->Solve(TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity);
    printf("%g %g %d %g\n", optimization->max_intracorrespondence_euclidean_distance, rigidity, 
      optimization->NCorrespondences(), optimization->Score());
    optimization->max_intracorrespondence_euclidean_distance *= 0.5;
    rigidity *= 0.5;
    break; }

    case GLUT_KEY_DOWN:
    case GLUT_KEY_UP:
    if (selected_image[side]) {
      GSVTapestry *tapestry = selected_image[side]->Tapestry();
      int index = selected_image[side]->TapestryIndex();
      if ((key == GLUT_KEY_DOWN) && (index > 0)) 
        selected_image[side] = tapestry->Image(index-1);
      else if ((key == GLUT_KEY_UP) && (index < tapestry->NImages()-1)) 
        selected_image[side] = tapestry->Image(index+1);
    }
    break;

    case GLUT_KEY_RIGHT:
    case GLUT_KEY_LEFT:
    if (selected_image[side]) {
      GSVPanorama *panorama = selected_image[side]->Panorama();
      int index = selected_image[side]->PanoramaIndex();
      if (key == GLUT_KEY_LEFT) index = (index > 0) ? index-1 : 7;
      if (key == GLUT_KEY_RIGHT) index = (index < 7) ? index+1 : 0;
      selected_image[side] = panorama->Image(index);
    }
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Remember side 
  GLUTside = side;

  // Redraw
  glutPostRedisplay();
}



////////////////////////////////////////////////////////////////////////
// Top-Level functions
////////////////////////////////////////////////////////////////////////

static void 
InitInterface(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA); 
  GLUTwindow = glutCreateWindow("GSV Street Viewer");

  // Initialize background color
  glClearColor(0, 0, 0, 0);

  // Initialize lights
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glEnable(GL_NORMALIZE);

  // Initialize headlight
  static GLfloat light0_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light0_specular[] = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat light0_position[] = { 0.0, 0.0, 1.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glEnable(GL_LIGHT0);

  // Initialize side light 
  static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light1_specular[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light1_position[] = { 0.5, -0.5, 0.866, 0 };
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
  glEnable(GL_LIGHT1);

  // Initialize materials 
  static GLfloat ambient_material[] = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat diffuse_material[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat specular_material[] = { 1, 1, 1, 1 };
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_material);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_material);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_material);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10);

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
  atexit(AtExit);
}



static void 
RunInterface(void)
{
  // Initialize center points
  center_point[0] = scene->BBox().Centroid();
  center_point[1] = scene->BBox().Centroid();

  // Initialize viewers
  for (int side = 0; side < 2; side++) {
    R3Box bbox = scene->BBox();
    assert(!bbox.IsEmpty());
    RNLength r = bbox.DiagonalRadius();
    assert((r > 0.0) && RNIsFinite(r));
    R3Point origin = center_point[side] + R3posz_vector * (2.5 * r);
    RNAngle yfov = RN_PI_OVER_FOUR;
    RNScalar window_aspect = (double) GLUTwindow_height / (double) GLUTwindow_width;
    RNAngle xfov = atan(tan(yfov) / window_aspect);
    R3Camera camera(origin, R3Vector(0, 0, -1), R3Vector(0, 1, 0), xfov, yfov, 0.001 * r, 1000.0 * r);
    R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
    viewer[side] = new R3Viewer(camera, viewport);
  }

  // Run main loop -- never returns
  glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-debug")) { print_debug = 1; }
      else if (!strcmp(*argv, "-interactive")) { interactive = 1; }
      else if (!strcmp(*argv, "-create_scan_pole_features")) { create_features = create_scan_pole_features = 1; }
      else if (!strcmp(*argv, "-create_scan_curb_features")) { create_features = create_scan_curb_features = 1; }
      else if (!strcmp(*argv, "-create_scan_edge_features")) { create_features = create_scan_edge_features = 1; }      
      else if (!strcmp(*argv, "-create_scan_plane_features")) { create_features = create_scan_plane_features = 1; }
      else if (!strcmp(*argv, "-create_image_corner_features")) { create_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_image_sift_features")) { create_features = create_image_sift_features = 1; }
      else if (!strcmp(*argv, "-create_image_line_features")) { create_features = create_image_line_features = 1; }
      else if (!strcmp(*argv, "-create_scan_features")) { create_features = create_scan_pole_features = create_scan_plane_features = 1; }
      else if (!strcmp(*argv, "-create_image_features")) { create_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_features")) { create_features = create_scan_pole_features = create_scan_plane_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_point_point_correspondences")) { create_correspondences = create_point_point_correspondences = 1; }
      else if (!strcmp(*argv, "-create_plane_plane_correspondences")) { create_correspondences = create_plane_plane_correspondences = 1; }
      else if (!strcmp(*argv, "-create_pixel_point_correspondences")) { create_correspondences = create_pixel_point_correspondences = 1; }
      else if (!strcmp(*argv, "-create_pixel_plane_correspondences")) { create_correspondences = create_pixel_plane_correspondences = 1; }
      else if (!strcmp(*argv, "-create_pixel_pixel_correspondences")) { create_correspondences = create_pixel_pixel_correspondences = 1; }
      else if (!strcmp(*argv, "-create_scan_correspondences")) { create_correspondences = create_point_point_correspondences = create_plane_plane_correspondences = 1; }
      else if (!strcmp(*argv, "-create_image_correspondences")) { create_correspondences = create_pixel_pixel_correspondences = 1; }
      else if (!strcmp(*argv, "-create_correspondences")) { create_correspondences = 
        create_point_point_correspondences = create_plane_plane_correspondences = create_pixel_pixel_correspondences = 1; }
        else if (!strcmp(*argv, "-create_transformations")) { create_transformations = 1; }
        else if (!strcmp(*argv, "-load_points_dynamically")) { load_points_dynamically = 1; }
        else if (!strcmp(*argv, "-load_points_statically")) { load_points_dynamically = 0; }
        else if (!strcmp(*argv, "-load_points_from_images")) { load_points_from_images = 1; }
        else if (!strcmp(*argv, "-update_pixel_features")) { update_pixel_features = 1; }
        else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene_name = *argv; }
        else if (!strcmp(*argv, "-output_scene")) { argc--; argv++; output_scene_name = *argv; }
        else if (!strcmp(*argv, "-input_features")) { argc--; argv++; input_correspondences_name = *argv; }
        else if (!strcmp(*argv, "-output_features")) { argc--; argv++; output_correspondences_name = *argv; }
        else if (!strcmp(*argv, "-input_correspondences")) { argc--; argv++; input_correspondences_name = *argv; }
        else if (!strcmp(*argv, "-output_correspondences")) { argc--; argv++; output_correspondences_name = *argv; }
        else if (!strcmp(*argv, "-input_transformations")) { argc--; argv++; input_transformations_name = *argv; }
        else if (!strcmp(*argv, "-output_transformations")) { argc--; argv++; output_transformations_name = *argv; }
        else if (!strcmp(*argv, "-output_error_plot")) { argc--; argv++; output_error_plot_name = *argv; }
        else if (!strcmp(*argv, "-rigidity")) { argc--; argv++; rigidity = atof(*argv); }
        else if (!strcmp(*argv, "-run_number")) { 
         argc--; argv++; 
         show_all_runs = 0; 
         run_number = atoi(*argv);
         if (run_number >= 2 && run_number <= 9) run_number += 10;
         else if (run_number >= 10 && run_number <= 19) run_number -= 8;
       }
       else if (!strcmp(*argv, "-segment_number")) { argc--; argv++; show_all_segments = 0; segment_number = atoi(*argv); }
       else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
       argv++; argc--;
     }
     else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_scene_name) output_scene_name = *argv;
      else if (!input_correspondences_name) input_correspondences_name = *argv;
      else if (!output_correspondences_name) output_correspondences_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsvalign input_scenefile [output_scenefile] [-interactive] [options]\n");
    return FALSE;
  }

  // Update point loading parameters
  if (!strstr(input_scene_name, ".gsv")) {
    load_points_dynamically = FALSE;
  }
  if (load_points_from_images) {
    load_points_dynamically = TRUE;
  }

  // Make interactive if no output
  if (!output_scene_name && !output_correspondences_name && !output_transformations_name && !output_error_plot_name) {
    interactive = TRUE;
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
  scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Allocate optimization data structures
  optimization = new GSVPoseOptimization(scene);
  if (!optimization) exit(-1);

  // Read correspondences
  if (input_correspondences_name) {
    if (!ReadCorrespondences(optimization, input_correspondences_name)) exit(-1);
  }

  // Read transformations
  if (input_transformations_name) {
    if (!ReadTransformations(optimization, input_transformations_name)) exit(-1);
  }

  // Create features
  if (create_features) {
    if (!CreateFeatures(optimization)) exit(-1);
  }

  // Update feature positions based on ray intersections
  if (update_pixel_features) {
    R3Box cull_box = R3null_box;
    if (!strcmp(input_scene_name, "test.gsv")) cull_box = R3Box(-RN_INFINITY, -600, -RN_INFINITY, -500, RN_INFINITY, RN_INFINITY);
    if (!UpdatePixelFeatures(optimization, cull_box)) exit(-1);
  }

  // Create correspondences
  if (create_correspondences) {
    if (!CreateCorrespondences(optimization)) exit(-1);
  }

  // Create transformations
  if (create_transformations) {
    if (!Solve(optimization, 4, solver)) exit(-1);
  }

  // Check if interactive
  if (interactive) {
    // Initialize GLUT interface
    InitInterface(&argc, argv);
    
    // Run GLUT interface
    RunInterface();
  }
  else {
    // Write correspondences
    if (output_correspondences_name) {
      if (!WriteCorrespondences(optimization, output_correspondences_name)) exit(-1);
    }

    // Write transformations
    if (output_transformations_name) {
      if (!WriteTransformations(optimization, output_transformations_name)) exit(-1);
    }

    // Write error plot
    if (output_error_plot_name) {
      if (!WriteErrorPlot(optimization, output_error_plot_name)) exit(-1);
    }

    // Write scene
    if (output_scene_name) {
      if (!WriteScene(scene, optimization, output_scene_name)) exit(-1);
    }
  }

  // Return success 
  return 0;
}


