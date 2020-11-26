// Source file for the google viewer program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static const char *input_scene_name = NULL;
static const char *input_image_directory = NULL;
static const char *input_initial_segmentation_basename = NULL;
static const char *input_seed_mask_basename = NULL;
static const char *input_inclusion_mask_basename = NULL;
static const char *output_image_directory = NULL;
static const char *output_image_basename = "Plane";
static int segment_SA_images = 0;
static int segment_DA_images = 0;
static int segment_DH_images = 0;
static int write_planar_grids = 0;
static int write_xy_spans = 0;
static int print_verbose = 0;
static int print_debug = 0;

// Masking parameters
static int mask_shadow_boundaries_from_inclusion = 0;
static int mask_boundaries_from_seeds = 1;
static int mask_inconsistent_normals_from_seeds = 1;
static int mask_sensor_vehicle = 1;

// Point pruning parameters
static int pca_index_radius = 2;
static double pca_min_lambda1 = 1E-2;
static double pca_min_lambda2 = 0;
static double pca_max_lambda2 = 0;
static double pca_max_lambda3 = 1E-2;
static double pca_max_lambda21 = 0;
static double pca_max_lambda31 = 0.25;
static double pca_max_lambda32 = 0.25;

// Hierarchical clustering parameters
static int max_neighbor_scanline_index_difference = 64;
static int max_neighbor_point_index_difference = 8;
static double max_neighbor_centroid_distance = 4;
static double max_neighbor_offplane_distance[3] = { 0.0, 0.0, 1.0 } ;
static double max_neighbor_axis_angle[3] = { 0.0, 0.0, RN_PI / 6.0 };
static double min_neighbor_affinity = 0;
static double max_pair_centroid_distance = 100;
static double max_pair_offplane_distance[3] = { 0.0, 0.0, 1.0 };
static double max_pair_offplane_angle[3] = { 0.0, 0.0, RN_PI / 8.0 };
static double max_pair_axis_angle[3] = { 0.0, 0.0, RN_PI / 8.0 };
static double min_pair_affinity = FLT_MIN;
static double min_pair_affinity_ratio = 0;



////////////////////////////////////////////////////////////////////////
// I/O Functions
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
  if (!scene->ReadFile(filename, FALSE)) {
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
// Planar Grid Functions
////////////////////////////////////////////////////////////////////////

static int
WritePlanarGrids(GSVScan *scan, 
  const char *input_image_directory, 
  const char *output_image_directory, const char *output_image_basename, 
  const char *image_parameterization, 
  int min_points = 16384, RNLength planar_grid_spacing = 0.1)
{
  // Get convenient variables
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Read position grids (from input image directory)
  char image_name[4096];
  R2Grid position_x_grid;
  R2Grid position_y_grid;
  R2Grid position_z_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_z_grid.Read(image_name)) return 0;

  // Read plane segment grids (from output image directory)
  R2Grid plane_segment_a_grid;
  R2Grid plane_segment_b_grid;
  R2Grid plane_segment_c_grid;
  R2Grid plane_segment_d_grid;
  R2Grid plane_segment_id_grid;
  R2Grid plane_segment_size_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneA.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_a_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneB.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_b_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneC.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_c_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneD.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_d_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneId.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_id_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_LargePlaneSize.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!plane_segment_size_grid.Read(image_name)) return 0;

  // Read color grids
  R2Grid color_red_grid, color_green_grid, color_blue_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorRed.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!color_red_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorGreen.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!color_green_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorBlue.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!color_blue_grid.Read(image_name)) return 0;

  // Read other attribute grids
  R2Grid height_grid, density_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Height.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!height_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Density.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!density_grid.Read(image_name)) return 0;

  // Parse segments
  RNArray<RNArray<const RNScalar *> *> segments;
  const RNScalar *grid_values = plane_segment_id_grid.GridValues();
  for (int i = 0; i < plane_segment_id_grid.NEntries(); i++) {
    int id_value = plane_segment_id_grid.GridValue(i);
    if (id_value == R2_GRID_UNKNOWN_VALUE) continue;
    int id = (int) (id_value + 0.5);
    while (id >= segments.NEntries()) segments.Insert(new RNArray<const RNScalar *>());
    RNArray<const RNScalar *> *segment = segments.Kth(id);
    segment->Insert(&grid_values[i]);
  }

  // Process segments
  for (int id = 0; id < segments.NEntries(); id++) {
    RNArray<const RNScalar *> *segment = segments.Kth(id);
    if (segment->NEntries() < min_points) continue;

    // Construct plane
    const RNScalar *grid_valuep = segment->Kth(0);
    int grid_index = grid_valuep - grid_values;
    RNCoord a = plane_segment_a_grid.GridValue(grid_index);
    if (a == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord b = plane_segment_b_grid.GridValue(grid_index);
    if (b == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord c = plane_segment_c_grid.GridValue(grid_index);
    if (c == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord d = plane_segment_d_grid.GridValue(grid_index);
    if (d == R2_GRID_UNKNOWN_VALUE) continue;
    R3Vector normal(a, b, c);
    normal.Normalize();
    R3Plane plane(normal, d);

    // Construct bounding box
    R3Box bbox = R3null_box;
    for (int j = 0; j < segment->NEntries(); j++) {
      const RNScalar *grid_valuep = segment->Kth(j);
      int grid_index = grid_valuep - grid_values;
      RNCoord x = position_x_grid.GridValue(grid_index);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord y = position_y_grid.GridValue(grid_index);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord z = position_z_grid.GridValue(grid_index);
      if (z == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point position(x, y, z);
      bbox.Union(position);
    }

    // Construct planar grids
    R3Point centroid = bbox.Centroid();
    R3PlanarGrid index_planar_grid(plane, bbox, centroid, R3posz_vector, planar_grid_spacing);
    R3PlanarGrid displacement_planar_grid(index_planar_grid);
    R3PlanarGrid overlap_planar_grid(index_planar_grid);
    R3PlanarGrid color_red_planar_grid(index_planar_grid);
    R3PlanarGrid color_green_planar_grid(index_planar_grid);
    R3PlanarGrid color_blue_planar_grid(index_planar_grid);
    R3PlanarGrid height_planar_grid(index_planar_grid);
    R3PlanarGrid density_planar_grid(index_planar_grid);
    for (int j = 0; j < segment->NEntries(); j++) {
      // Get grid position
      int ix, iy;
      const RNScalar *grid_valuep = segment->Kth(j);
      int grid_index = grid_valuep - grid_values;
      plane_segment_id_grid.IndexToIndices(grid_index, ix, iy);
      if (ix >= plane_segment_id_grid.XResolution()-1) continue;
      if (iy >= plane_segment_id_grid.YResolution()-1) continue;

      // Get world position
      R3Point p00(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
      RNScalar id00 = plane_segment_id_grid.GridValue(ix, iy);
      assert((int) (id00 + 0.5) == id);
      if (id00 == R2_GRID_UNKNOWN_VALUE) continue;
      p00[0] = position_x_grid.GridValue(ix, iy);
      p00[1] = position_y_grid.GridValue(ix, iy);
      p00[2] = position_z_grid.GridValue(ix, iy);
      if (p00[0] == R2_GRID_UNKNOWN_VALUE) continue;
      if (p00[1] == R2_GRID_UNKNOWN_VALUE) continue;
      if (p00[2] == R2_GRID_UNKNOWN_VALUE) continue;

      // Get point on right
      R3Point p10(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
      RNScalar id10 = plane_segment_id_grid.GridValue(ix+1, iy);
      if ((int) (id10 + 0.5) == id) {
        p10[0] = position_x_grid.GridValue(ix+1, iy);
        p10[1] = position_y_grid.GridValue(ix+1, iy);
        p10[2] = position_z_grid.GridValue(ix+1, iy);
        if (p10[1] == R2_GRID_UNKNOWN_VALUE) p10[0] = R2_GRID_UNKNOWN_VALUE;
        if (p10[2] == R2_GRID_UNKNOWN_VALUE) p10[0] = R2_GRID_UNKNOWN_VALUE;
      }
      
      // Get point above 
      R3Point p01(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
      RNScalar id01 = plane_segment_id_grid.GridValue(ix, iy+1);
      if ((int) (id01 + 0.5) == id) {
        p01[0] = position_x_grid.GridValue(ix, iy+1);
        p01[1] = position_y_grid.GridValue(ix, iy+1);
        p01[2] = position_z_grid.GridValue(ix, iy+1);
        if (p01[1] == R2_GRID_UNKNOWN_VALUE) p01[0] = R2_GRID_UNKNOWN_VALUE;
        if (p01[2] == R2_GRID_UNKNOWN_VALUE) p01[0] = R2_GRID_UNKNOWN_VALUE;
      }
      
      // Get point diagonally right and above
      R3Point p11(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
      RNScalar id11 = plane_segment_id_grid.GridValue(ix+1, iy+1);
      if ((int) (id11 + 0.5) == id) {
        p11[0] = position_x_grid.GridValue(ix+1, iy+1);
        p11[1] = position_y_grid.GridValue(ix+1, iy+1);
        p11[2] = position_z_grid.GridValue(ix+1, iy+1);
        if (p11[1] == R2_GRID_UNKNOWN_VALUE) p11[0] = R2_GRID_UNKNOWN_VALUE;
        if (p11[2] == R2_GRID_UNKNOWN_VALUE) p11[0] = R2_GRID_UNKNOWN_VALUE;
      }

      // Draw triangles
      if ((p00[0] != R2_GRID_UNKNOWN_VALUE) && (p11[0] != R2_GRID_UNKNOWN_VALUE)) {
        if (p10[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d00 = R3SignedDistance(plane, p00);
          RNScalar d10 = R3SignedDistance(plane, p10);
          RNScalar d11 = R3SignedDistance(plane, p11);
          RNScalar r00 = color_red_grid.GridValue(ix, iy);
          RNScalar r10 = color_red_grid.GridValue(ix+1, iy);
          RNScalar r11 = color_red_grid.GridValue(ix+1, iy+1);
          RNScalar g00 = color_green_grid.GridValue(ix, iy);
          RNScalar g10 = color_green_grid.GridValue(ix+1, iy);
          RNScalar g11 = color_green_grid.GridValue(ix+1, iy+1);
          RNScalar b00 = color_blue_grid.GridValue(ix, iy);
          RNScalar b10 = color_blue_grid.GridValue(ix+1, iy);
          RNScalar b11 = color_blue_grid.GridValue(ix+1, iy+1);
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p10, p11, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p10, p11, d00, d10, d11);
          color_red_planar_grid.RasterizeWorldTriangle(p00, p10, p11, r00, r10, r11);
          color_green_planar_grid.RasterizeWorldTriangle(p00, p10, p11, g00, g10, g11);
          color_blue_planar_grid.RasterizeWorldTriangle(p00, p10, p11, b00, b10, b11);
          height_planar_grid.RasterizeWorldTriangle(p00, p10, p11, h00, h10, h11);
          density_planar_grid.RasterizeWorldTriangle(p00, p10, p11, s00, s10, s11);
          overlap_planar_grid.RasterizeWorldTriangle(p00, p10, p11, 1.0);
        }
        if (p01[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d00 = R3SignedDistance(plane, p00);
          RNScalar d01 = R3SignedDistance(plane, p01);
          RNScalar d11 = R3SignedDistance(plane, p11);
          RNScalar r00 = color_red_grid.GridValue(ix, iy);
          RNScalar r01 = color_red_grid.GridValue(ix, iy+1);
          RNScalar r11 = color_red_grid.GridValue(ix+1, iy+1);
          RNScalar g00 = color_green_grid.GridValue(ix, iy);
          RNScalar g01 = color_green_grid.GridValue(ix, iy+1);
          RNScalar g11 = color_green_grid.GridValue(ix+1, iy+1);
          RNScalar b00 = color_blue_grid.GridValue(ix, iy);
          RNScalar b01 = color_blue_grid.GridValue(ix, iy+1);
          RNScalar b11 = color_blue_grid.GridValue(ix+1, iy+1);
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p11, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p11, p01, d00, d11, d01);
          color_red_planar_grid.RasterizeWorldTriangle(p00, p11, p01, r00, r11, r01);
          color_green_planar_grid.RasterizeWorldTriangle(p00, p11, p01, g00, g11, g01);
          color_blue_planar_grid.RasterizeWorldTriangle(p00, p11, p01, b00, b11, b01);
          height_planar_grid.RasterizeWorldTriangle(p00, p11, p01, h00, h11, h01);
          density_planar_grid.RasterizeWorldTriangle(p00, p11, p01, s00, s11, s01);
          overlap_planar_grid.RasterizeWorldTriangle(p00, p11, p01, 1.0);
        }
      }
      else if ((p10[0] != R2_GRID_UNKNOWN_VALUE) && (p01[0] != R2_GRID_UNKNOWN_VALUE)) {
        if (p00[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d00 = R3SignedDistance(plane, p00);
          RNScalar d01 = R3SignedDistance(plane, p01);
          RNScalar d10 = R3SignedDistance(plane, p10);
          RNScalar r00 = color_red_grid.GridValue(ix, iy);
          RNScalar r10 = color_red_grid.GridValue(ix+1, iy);
          RNScalar r01 = color_red_grid.GridValue(ix, iy+1);
          RNScalar g00 = color_green_grid.GridValue(ix, iy);
          RNScalar g10 = color_green_grid.GridValue(ix+1, iy);
          RNScalar g01 = color_green_grid.GridValue(ix, iy+1);
          RNScalar b00 = color_blue_grid.GridValue(ix, iy);
          RNScalar b10 = color_blue_grid.GridValue(ix+1, iy);
          RNScalar b01 = color_blue_grid.GridValue(ix, iy+1);
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p10, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p10, p01, d00, d10, d01);
          color_red_planar_grid.RasterizeWorldTriangle(p00, p10, p01, r00, r10, r01);
          color_green_planar_grid.RasterizeWorldTriangle(p00, p10, p01, g00, g10, g01);
          color_blue_planar_grid.RasterizeWorldTriangle(p00, p10, p01, b00, b10, b01);
          height_planar_grid.RasterizeWorldTriangle(p00, p10, p01, h00, h10, h01);
          density_planar_grid.RasterizeWorldTriangle(p00, p10, p01, s00, s10, s01);
          overlap_planar_grid.RasterizeWorldTriangle(p00, p10, p01, 1.0);
        }
        if (p11[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d01 = R3SignedDistance(plane, p01);
          RNScalar d10 = R3SignedDistance(plane, p10);
          RNScalar d11 = R3SignedDistance(plane, p11);
          RNScalar r10 = color_red_grid.GridValue(ix+1, iy);
          RNScalar r01 = color_red_grid.GridValue(ix, iy+1);
          RNScalar r11 = color_red_grid.GridValue(ix+1, iy+1);
          RNScalar g10 = color_green_grid.GridValue(ix+1, iy);
          RNScalar g01 = color_green_grid.GridValue(ix, iy+1);
          RNScalar g11 = color_green_grid.GridValue(ix+1, iy+1);
          RNScalar b10 = color_blue_grid.GridValue(ix+1, iy);
          RNScalar b01 = color_blue_grid.GridValue(ix, iy+1);
          RNScalar b11 = color_blue_grid.GridValue(ix+1, iy+1);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p10, p11, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p10, p11, p01, d10, d11, d01);
          color_red_planar_grid.RasterizeWorldTriangle(p10, p11, p01, r10, r11, r01);
          color_green_planar_grid.RasterizeWorldTriangle(p10, p11, p01, g10, g11, g01);
          color_blue_planar_grid.RasterizeWorldTriangle(p10, p11, p01, b10, b11, b01);
          height_planar_grid.RasterizeWorldTriangle(p10, p11, p01, h10, h11, h01);
          density_planar_grid.RasterizeWorldTriangle(p10, p11, p01, s10, s11, s01);
          overlap_planar_grid.RasterizeWorldTriangle(p10, p11, p01, 1.0);
        }
      }
    }

    // Mark unknown areas
    R3PlanarGrid mask(overlap_planar_grid);
    mask.Threshold(0.5, R2_GRID_UNKNOWN_VALUE, 1);
    displacement_planar_grid.Mask(mask);
    color_red_planar_grid.Mask(mask);
    color_green_planar_grid.Mask(mask);
    color_blue_planar_grid.Mask(mask);
    height_planar_grid.Mask(mask);
    density_planar_grid.Mask(mask);
    index_planar_grid.Mask(mask);

    // Compute means
    displacement_planar_grid.Divide(overlap_planar_grid);
    color_red_planar_grid.Divide(overlap_planar_grid);
    color_green_planar_grid.Divide(overlap_planar_grid);
    color_blue_planar_grid.Divide(overlap_planar_grid);
    height_planar_grid.Divide(overlap_planar_grid);
    density_planar_grid.Divide(overlap_planar_grid);
    index_planar_grid.Divide(overlap_planar_grid); // NOT RIGHT

    // Create color image
    R2Image rgb_image(color_red_planar_grid.XResolution(), color_red_planar_grid.YResolution());
    for (int ix = 0; ix < color_red_planar_grid.XResolution(); ix++) {
      for (int iy = 0; iy < color_red_planar_grid.YResolution(); iy++) {
        RNScalar r = color_red_planar_grid.GridValue(ix, iy);
        if (r == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar g = color_green_planar_grid.GridValue(ix, iy);
        if (g == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar b = color_blue_planar_grid.GridValue(ix, iy);
        if (b == R2_GRID_UNKNOWN_VALUE) continue;
        RNRgb color(r, g, b);
        rgb_image.SetPixelRGB(ix, iy, color);
      }
    }

    // Write planar grids
    char planar_grid_name[1024];
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_Displacement.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!displacement_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_%sGridIndex.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id, image_parameterization);
    if (!index_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_ColorRed.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!color_red_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_ColorGreen.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!color_green_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_ColorBlue.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!color_blue_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_Color.bmp", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    rgb_image.Write(planar_grid_name);
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_Height.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!height_planar_grid.WriteFile(planar_grid_name)) return 0;
    sprintf(planar_grid_name, "%s/%s/%02d_%02d_%06d_ST_Density.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!density_planar_grid.WriteFile(planar_grid_name)) return 0;

    // Print debug message
    if (print_debug) {
      printf("  %20s %2d %2d : %6.3f %6.3f %6.3f %6.3f : %6.3f %6.3f %6.3f : %d\n", 
        run->Name(), segment_index, scan_index,  
        plane.A(), plane.B(), plane.C(), plane.D(), 
        bbox.XLength(), bbox.YLength(), bbox.ZLength(),
        segment->NEntries());
    }
  }

  // Delete segments
  for (int i = 0; i < segments.NEntries(); i++) {
    RNArray<const RNScalar *> *segment = segments.Kth(i);
    delete segment;
  }

  // Return success
  return 1;
}



static int
WritePlanarGrids(GSVScene *scene, 
  const char *input_image_directory, 
  const char *output_image_directory, const char *output_image_basename, 
  const char *image_parameterization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Writing %s planar grids ...\n", image_parameterization);
    fflush(stdout);
  }

  // Write planar grids for all scans
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (!WritePlanarGrids(scan, input_image_directory, 
          output_image_directory, output_image_basename, 
          image_parameterization)) return 0;
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
// Span Functions
////////////////////////////////////////////////////////////////////////

static int
WriteXYSpans(FILE *fp, GSVScan *scan, 
  const char *input_image_directory, const char *output_image_directory, 
  const char *output_image_basename, const char *image_parameterization, 
  int min_points = 16, RNScalar max_normal_z = 0.1)
{
  // Get convenient variables
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Read position grids (from input image directory)
  char image_name[4096];
  R2Grid position_x_grid;
  R2Grid position_y_grid;
  R2Grid position_z_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_z_grid.Read(image_name)) return 0;

  // Read plane segment grids (from output image directory)
  R2Grid plane_segment_a_grid;
  R2Grid plane_segment_b_grid;
  R2Grid plane_segment_c_grid;
  R2Grid plane_segment_d_grid;
  R2Grid plane_segment_id_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sA.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, output_image_basename);
  if (!plane_segment_a_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sB.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, output_image_basename);
  if (!plane_segment_b_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sC.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, output_image_basename);
  if (!plane_segment_c_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sD.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, output_image_basename);
  if (!plane_segment_d_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sId.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization, output_image_basename);
  if (!plane_segment_id_grid.Read(image_name)) return 0;

  // Parse segments
  RNArray<RNArray<const RNScalar *> *> segments;
  const RNScalar *grid_values = plane_segment_id_grid.GridValues();
  for (int i = 0; i < plane_segment_id_grid.NEntries(); i++) {
    int id_value = plane_segment_id_grid.GridValue(i);
    if (id_value == R2_GRID_UNKNOWN_VALUE) continue;
    int id = (int) (id_value + 0.5);
    while (id >= segments.NEntries()) segments.Insert(new RNArray<const RNScalar *>());
    RNArray<const RNScalar *> *segment = segments.Kth(id);
    segment->Insert(&grid_values[i]);
  }

  // Process segments
  for (int id = 0; id < segments.NEntries(); id++) {
    RNArray<const RNScalar *> *segment = segments.Kth(id);

    // Check number of points
    if (segment->NEntries() < min_points) continue;

    // Get plane info
    const RNScalar *grid_valuep = segment->Kth(0);
    int grid_index = grid_valuep - grid_values;
    RNCoord a = plane_segment_a_grid.GridValue(grid_index);
    if (a == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord b = plane_segment_b_grid.GridValue(grid_index);
    if (b == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord c = plane_segment_c_grid.GridValue(grid_index);
    if (c == R2_GRID_UNKNOWN_VALUE) continue;
    RNCoord d = plane_segment_d_grid.GridValue(grid_index);
    if (d == R2_GRID_UNKNOWN_VALUE) continue;

    // Check plane normal
    if (RNIsNotZero(c, max_normal_z)) continue;

    // Analyze points
    int count = 0;
    R2Point centroid(0, 0);
    for (int j = 0; j < segment->NEntries(); j++) {
      const RNScalar *grid_valuep = segment->Kth(j);
      int grid_index = grid_valuep - grid_values;
      RNCoord x = position_x_grid.GridValue(grid_index);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord y = position_y_grid.GridValue(grid_index);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord z = position_z_grid.GridValue(grid_index);
      if (z == R2_GRID_UNKNOWN_VALUE) continue;
      R2Point position(x, y);
      centroid += position;
      count++;
    }

    // Compute centroid
    if (count == 0) continue;
    if (count < min_points) continue;
    centroid /= count;

    // Compute 2D line on XY plane
    R2Vector normal(a, b);
    R2Vector direction(-b, a);
    direction.Normalize();
    R2Point start(centroid.X(), centroid.Y());
    R2Ray ray(start, direction);

    // Compute extent
    RNScalar max_d = 0;
    RNScalar min_t = FLT_MAX;
    RNScalar max_t = -FLT_MAX;
    RNScalar min_z = FLT_MAX;
    RNScalar max_z = -FLT_MAX;
    for (int j = 0; j < segment->NEntries(); j++) {
      const RNScalar *grid_valuep = segment->Kth(j);
      int grid_index = grid_valuep - grid_values;
      RNCoord x = position_x_grid.GridValue(grid_index);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord y = position_y_grid.GridValue(grid_index);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord z = position_z_grid.GridValue(grid_index);
      if (z == R2_GRID_UNKNOWN_VALUE) continue;
      R2Point p(x, y);
      RNScalar d = R2Distance(p, ray.Line());
      RNScalar t = ray.T(p);
      if (z < min_z) min_z = z;
      if (z > max_z) max_z = z;
      if (t < min_t) min_t = t;
      if (t > max_t) max_t = t;
      if (d > max_d) max_d = d;
    }

    // Write span
    R2Point p1 = ray.Point(min_t);
    R2Point p2 = ray.Point(max_t);
    RNLength length = R2Distance(p1, p2);
    if (RNIsZero(length)) continue;
    RNLength height = max_z - min_z;
    if (RNIsZero(height)) continue;
    RNArea area = length * height;
    RNScalar density = count / area;
    fprintf(fp, "  %20s %2d %2d : %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f   %9.3f %9.3f %9.3f %9.3f %9d %9.3f\n", 
      run->Name(), segment_index, scan_index,  p1.X(), p1.Y(), p2.X(), p2.Y(), min_z, max_z, 
      max_d, length, height, area, count, density);
  }

  // Delete segments
  for (int i = 0; i < segments.NEntries(); i++) {
    RNArray<const RNScalar *> *segment = segments.Kth(i);
    delete segment;
  }

  // Return success
  return 1;
}



static int
WriteXYSpans(GSVScene *scene, 
  const char *input_image_directory, 
  const char *output_image_directory, const char *output_image_basename, 
  const char *image_parameterization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Writing %s spans ...\n", image_parameterization);
    fflush(stdout);
  }

  // Create filename
  char filename[4096];
  sprintf(filename, "%s/XY_%s_%sSpans.txt", output_image_directory, image_parameterization, output_image_basename);
  if (RNFileExists(filename)) return 1;

  // Open output file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open span file %s\n", filename);
    return 0;
  }

  // Write planar grids for all scans
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (!WriteXYSpans(fp, scan, input_image_directory, output_image_directory, 
          output_image_basename, image_parameterization)) return 0;
      }
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("  Done in %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Process normal images
////////////////////////////////////////////////////////////////////////

static void
MaskPointsOnSensorVehicle(R2Grid& mask)
{
  // Mask positions that are on sensor car 
  for (int i = 0; i < mask.XResolution(); i++) {
    for (int j = 0; j < 20; j++) {
      mask.SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
    }
  }
}



static void
MaskBoundaries(R2Grid& mask, const R2Grid& boundary_type_grid, int masked_boundary_types)
{
  // Mask boundaries of types indicated by bit mask
  for (int i = 0; i < mask.XResolution(); i++) {
    for (int j = 0; j < mask.YResolution(); j++) {
      if (mask.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
      if (mask.GridValue(i, j) == 0) continue;
      RNScalar boundary_type_value = boundary_type_grid.GridValue(i, j);
      int boundary_type = (int) (boundary_type_value + 0.5);
      if ((boundary_type_value == R2_GRID_UNKNOWN_VALUE) || (boundary_type & masked_boundary_types)) {
        mask.SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
      }
    }
  }
}



static void
MaskOcclusionBoundaries(R2Grid& mask, const R2Grid& viewpoint_depth_grid, 
  RNLength max_depth_difference = 1, RNScalar min_depth_difference_ratio = 0.5)
{
  // Mask points that are on boundaries of depth discontinuities
  for (int i = 1; i < mask.XResolution()-1; i++) {
    for (int j = 1; j < mask.YResolution()-1; j++) {
      // Get depth
      RNScalar depth = viewpoint_depth_grid.GridValue(i, j);
      if (depth == R2_GRID_UNKNOWN_VALUE) continue;

      // Consider opposite pairs of neighbors
      for (int k = 1; k <= 3; k++) {
        int dix = k / 2;
        int diy = k % 2;
        RNScalar depthA = viewpoint_depth_grid.GridValue(i-dix, j-diy);
        if (depthA == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar depthB = viewpoint_depth_grid.GridValue(i+dix, j+diy);
        if (depthB == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar deltaA = fabs(depth - depthA);
        RNScalar deltaB = fabs(depth - depthB);
        RNScalar minAB = (deltaA > deltaB) ? deltaB : deltaA;
        RNScalar maxAB = (deltaA > deltaB) ? deltaA : deltaB;
        if (maxAB > max_depth_difference) {
          RNScalar ratio = minAB / maxAB;
          if (ratio < min_depth_difference_ratio) {
            mask.SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
            break;
          }
        }
      }
    }
  }
}



static void
MaskInconsistentNormals(R2Grid& mask,
  const R2Grid& normal_x_grid, const R2Grid& normal_y_grid, const R2Grid& normal_z_grid, 
  RNAngle angle_threshold = RN_PI / 8.0, RNScalar consistent_fraction_threshold = 0.85, int neighborhood_radius = 1)
{
  // Get convenient variables
  int xres = normal_x_grid.XResolution();
  int yres = normal_x_grid.YResolution();
  if (neighborhood_radius > xres-1) neighborhood_radius = xres-1;
  if (neighborhood_radius > yres-1) neighborhood_radius = yres-1;
  int neighborhood_diameter = 2*neighborhood_radius + 1;
  int count_threshold = (int) (consistent_fraction_threshold * (neighborhood_diameter * neighborhood_diameter - 1));

  // Create copies of grids
  R2Grid original_normal_x_grid(normal_x_grid);
  R2Grid original_normal_y_grid(normal_y_grid);
  R2Grid original_normal_z_grid(normal_z_grid);

  // Mask normals on top and bottom boundaries
  for (int i = 0; i < normal_x_grid.XResolution(); i++) { 
    for (int k = 0; k <= neighborhood_radius; k++) {
      mask.SetGridValue(i, k, R2_GRID_UNKNOWN_VALUE);
      mask.SetGridValue(i, yres-1-k, R2_GRID_UNKNOWN_VALUE);
    }
  }

  // Mask normals on left and right boundaries
  for (int j = 0; j < normal_x_grid.YResolution(); j++) { 
    for (int k = 0; k <= neighborhood_radius; k++) {
      mask.SetGridValue(k, j, R2_GRID_UNKNOWN_VALUE);
      mask.SetGridValue(xres-1-k, j, R2_GRID_UNKNOWN_VALUE);
    }
  }

  // Mask normals that are inconsistent with their neighbors
  for (int i = neighborhood_radius; i < normal_x_grid.XResolution()-neighborhood_radius; i++) {
    for (int j = neighborhood_radius; j < normal_x_grid.YResolution()-neighborhood_radius; j++) {
      RNScalar nx0 = original_normal_x_grid.GridValue(i, j);
      RNScalar ny0 = original_normal_y_grid.GridValue(i, j);
      RNScalar nz0 = original_normal_z_grid.GridValue(i, j);
      if (nx0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (ny0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (nz0 == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector n0(nx0, ny0, nz0);

      // Count neighbors with similar normals
      int count = 0;
      for (int s = -neighborhood_radius; s <= neighborhood_radius; s++) {
        for (int t = -neighborhood_radius; t <= neighborhood_radius; t++) {
          if ((s == 0) && (t == 0)) continue;

          // Get normal at adjacent point
          RNScalar nx1 = original_normal_x_grid.GridValue(i+s, j+t);
          RNScalar ny1 = original_normal_y_grid.GridValue(i+s, j+t);
          RNScalar nz1 = original_normal_z_grid.GridValue(i+s, j+t);
          if (nx1 == R2_GRID_UNKNOWN_VALUE) continue;
          if (ny1 == R2_GRID_UNKNOWN_VALUE) continue;
          if (nz1 == R2_GRID_UNKNOWN_VALUE) continue;
          R3Vector n1(nx1, ny1, nz1);

          // Compute angle between normals
          RNScalar dot = n1.Dot(n0);
          if (dot < 0) dot = -dot;
          RNAngle angle = (dot < 1) ? acos(dot) : 0;

          // Increment counter if within threshold
          if (angle < angle_threshold) count++;
        }
      }

      // Check count of neighbors with similar normals
      if (count < count_threshold) {
        mask.SetGridValue(i, j, R2_GRID_UNKNOWN_VALUE);
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Hierarchical clustering
////////////////////////////////////////////////////////////////////////

struct Cluster {
  Cluster *parent;
  RNArray<struct ClusterPair *> pairs;
  RNArray<const RNScalar *> points;
  int initial_segment_membership;
  R3Point centroid;
  R3Vector axes[3];
  RNScalar variances[3];
  RNScalar affinity;
};

struct ClusterPair {
  Cluster *clusters[2];
  RNScalar affinity;
  ClusterPair **heapentry;
};


static int
CompareClusters(const void *data1, const void *data2)
{
  Cluster *cluster1 = *((Cluster **) data1);
  Cluster *cluster2 = *((Cluster **) data2);
  return (cluster2->points.NEntries() - cluster1->points.NEntries());
}


static RNScalar
NeighborAffinity(Cluster *cluster0, Cluster *cluster1)
{
  // Check initial segment memberships
  if (cluster0->initial_segment_membership && cluster1->initial_segment_membership && 
    !(cluster0->initial_segment_membership & cluster1->initial_segment_membership)) return - 1.0;

  // Compute centroid distance factor
  RNScalar centroid_distance_factor = 1;
  if (max_neighbor_centroid_distance > 0) {
    RNLength centroid_distance = R3Distance(cluster0->centroid, cluster1->centroid);
    if (centroid_distance > max_neighbor_centroid_distance) return -1;
    centroid_distance_factor = 1.0 - centroid_distance / max_neighbor_centroid_distance;
  }

  // Compute offplane distance factors
  RNScalar offplane0_distance_factor = 1;
  RNScalar offplane1_distance_factor = 1;
  for (int i = 0; i < 3; i++) {
    if (max_neighbor_offplane_distance[i] > 0) {
      // Compute point1-plane0 distance affinity
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane0_distance_factor *= 0.001;
      }
      else {
        R3Plane plane0(cluster0->centroid, cluster0->axes[i]);
        RNLength offplane0_distance = R3Distance(plane0, cluster1->centroid);
        if (offplane0_distance > max_neighbor_offplane_distance[i]) return -1;
        offplane0_distance_factor *= 1.0 - offplane0_distance / max_neighbor_offplane_distance[i];
      }
      
      // Compute point0-plane1 distance affinity
      if (cluster1->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane1_distance_factor *= 0.001;
      }
      else {
        R3Plane plane1(cluster1->centroid, cluster1->axes[i]);
        RNLength offplane1_distance = R3Distance(plane1, cluster0->centroid);
        if (offplane1_distance > max_neighbor_offplane_distance[i]) return -1;
        offplane1_distance_factor *= 1.0 - offplane1_distance / max_neighbor_offplane_distance[i];
      }
    }
  }

  // Compute axis angle factors
  RNScalar axis_angle_factor = 1;
  for (int i = 0; i < 3; i++) {
    if (max_neighbor_axis_angle[i] > 0) {
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        axis_angle_factor *= 0.001;
      }
      else if (cluster1->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        axis_angle_factor *= 0.001;
      }
      else {
        RNScalar dot = fabs(cluster0->axes[i].Dot(cluster1->axes[i]));
        RNAngle axis_angle = (dot < 1) ? acos(dot) : 0;
        if (axis_angle > max_neighbor_axis_angle[i]) return -1;
        axis_angle_factor *= 1.0 - axis_angle / max_neighbor_axis_angle[i];
      }
    }
  }

  // Compute affinity
  RNScalar affinity = 1;
  affinity *= centroid_distance_factor;
  affinity *= offplane0_distance_factor;
  affinity *= offplane1_distance_factor;
  affinity *= axis_angle_factor;

  // Return affinity
  return affinity;
}



static RNScalar
ClusterAffinity(Cluster *cluster0, Cluster *cluster1, RNBoolean print_debug = FALSE)
{
  // Check initial segment memberships
  if (cluster0->initial_segment_membership && cluster1->initial_segment_membership && 
    !(cluster0->initial_segment_membership & cluster1->initial_segment_membership)) return - 1.0;

  // Compute centroid distance factor
  RNScalar centroid_distance_factor = 1;
  if (max_pair_centroid_distance > 0) {
    RNLength centroid_distance = R3Distance(cluster0->centroid, cluster1->centroid);
    if (centroid_distance > max_pair_centroid_distance) return -1;
    centroid_distance_factor = 1.0 - centroid_distance / max_pair_centroid_distance;
    centroid_distance_factor = pow(centroid_distance_factor, 10);
  }

  // Compute offplane distance factors
  RNScalar offplane0_distance_factor = 1;
  RNScalar offplane1_distance_factor = 1;
  for (int i = 0; i < 3; i++) {
    if (max_pair_offplane_distance[i] > 0) {
      // Compute point1-plane0 distance
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane0_distance_factor *= 0.001;
      }
      else {
        R3Plane plane0(cluster0->centroid, cluster0->axes[i]);
        RNLength offplane0_distance = R3Distance(plane0, cluster1->centroid);
        if (offplane0_distance > max_pair_offplane_distance[i]) return -1;
        offplane0_distance_factor *= 1.0 - offplane0_distance / max_pair_offplane_distance[i];
      }

      // Compute point0-plane1 distance
      if (cluster1->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane1_distance_factor *= 0.001;
      }
      else {
        R3Plane plane1(cluster1->centroid, cluster1->axes[i]);
        RNLength offplane1_distance = R3Distance(plane1, cluster0->centroid);
        if (offplane1_distance > max_pair_offplane_distance[i]) return -1;
        offplane1_distance_factor *= 1.0 - offplane1_distance / max_pair_offplane_distance[i];
      }
    }
  }

  // Compute offplane angle factors
  RNScalar offplane0_angle_factor = 1;
  RNScalar offplane1_angle_factor = 1;
  for (int i = 0; i < 3; i++) {
    if (max_pair_offplane_angle[i] > 0) {
      // Compute point1-plane0 angle
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane0_angle_factor *= 0.001;
      }
      else {
        R3Vector v = cluster1->centroid - cluster0->centroid;
        RNLength d = v.Length();
        if (RNIsPositive(d) > 0) {
          RNScalar cos_angle = fabs(v.Dot(cluster0->axes[i])) / d;
          RNLength offplane0_angle = RN_PI_OVER_TWO - ((cos_angle < 1) ? acos(cos_angle) : 0);
          if (offplane0_angle > max_pair_offplane_angle[i]) return -1;
          offplane0_angle_factor *= 1.0 - offplane0_angle / max_pair_offplane_angle[i];
        }
      }

      // Compute point1-plane0 angle
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        offplane1_angle_factor *= 0.001;
      }
      else {
        R3Vector v = cluster0->centroid - cluster1->centroid;
        RNLength d = v.Length();
        if (RNIsPositive(d) > 0) {
          RNScalar cos_angle = fabs(v.Dot(cluster1->axes[i])) / d;
          RNLength offplane1_angle = RN_PI_OVER_TWO - ((cos_angle < 1) ? acos(cos_angle) : 0);
          if (offplane1_angle > max_pair_offplane_angle[i]) return -1;
          offplane1_angle_factor *= 1.0 - offplane1_angle / max_pair_offplane_angle[i];
        }
      }
    }
  }

  // Compute normal angle factor
  RNScalar axis_angle_factor = 1;
  for (int i = 0; i < 3; i++) {
    if (max_pair_axis_angle[i] > 0) {
      if (cluster0->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        axis_angle_factor *= 0.001;
      }
      else if (cluster1->axes[i][0] == R2_GRID_UNKNOWN_VALUE) {
        axis_angle_factor *= 0.001;
      }
      else {
        RNScalar dot = fabs(cluster0->axes[i].Dot(cluster1->axes[i]));
        RNAngle axis_angle = (dot < 1) ? acos(dot) : 0;
        if (axis_angle > max_pair_axis_angle[i]) return -1;
        axis_angle_factor *= 1.0 - axis_angle / max_pair_axis_angle[i];
      }
    }
  }

  // Compute cluster size factor (favor merges of small clusters into big ones)
  int points0 = cluster0->points.NEntries();
  int points1 = cluster1->points.NEntries();
  RNScalar max_points = (points0 > points1) ? points0 : points1;
  RNScalar total_points = points0 + points1;
  RNScalar cluster_size_factor = max_points / total_points;

  // Compute affinity
  RNScalar affinity = 1;
  affinity *= centroid_distance_factor;
  affinity *= offplane0_distance_factor;
  affinity *= offplane1_distance_factor;
  affinity *= offplane0_angle_factor;
  affinity *= offplane1_angle_factor;
  affinity *= axis_angle_factor;
  affinity *= cluster_size_factor;

  // Print debug statement
  if (print_debug) {
    printf("  A %6.3f : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", affinity,
      centroid_distance_factor, offplane0_distance_factor, offplane1_distance_factor,
      offplane0_angle_factor, offplane1_angle_factor, axis_angle_factor, cluster_size_factor);
  }

  // Return affinity
  return affinity;
}



static void
UpdateClusterTriad(Cluster *cluster, 
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  RNScalar min_variance = 0.01) 
{
  // Check number of points
  if (cluster->points.NEntries() < 5) return;

  // Create array of point positions
  int npositions = 0;
  const int max_positions = 1024;
  static R3Point positions[max_positions];
  int step = 1 + cluster->points.NEntries() / max_positions;
  for (int i = 0; i < cluster->points.NEntries(); i += step) {
    const RNScalar *point = cluster->points.Kth(i);
    int point_index = point - position_x_grid.GridValues();
    RNScalar x = position_x_grid.GridValue(point_index);
    RNScalar y = position_y_grid.GridValue(point_index);
    RNScalar z = position_z_grid.GridValue(point_index);
    positions[npositions++] = R3Point(x, y, z);
    if (npositions >= max_positions) break;
  }

  // Compute triad and variances from point positions
  RNScalar variances[3];
  R3Point centroid = R3Centroid(npositions, positions);
  R3Triad triad = R3PrincipleAxes(centroid, npositions, positions, NULL, variances);

  // Update axes[0]
  cluster->axes[1] = triad[0];

  // Update axes[1]
  if ((cluster->axes[1][0] == R2_GRID_UNKNOWN_VALUE) || (variances[0] > min_variance)) {
    cluster->axes[1] = triad[1];
  }

  // Update axes[2]
  if ((cluster->axes[2][0] == R2_GRID_UNKNOWN_VALUE) || (variances[1] > min_variance)) {
    cluster->axes[2] = triad[2];
  }

  // Update variances
  cluster->variances[0] = variances[0];
  cluster->variances[1] = variances[1];
  cluster->variances[2] = variances[2];
}



static int
InitializeCluster(GSVScan *scan, Cluster *cluster, int grid_index, 
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  const R2Grid& normal_x_grid, const R2Grid& normal_y_grid, const R2Grid& normal_z_grid,
  const R2Grid *initial_segmentation_grid = NULL, const R2Grid *seed_mask_grid = NULL, const R2Grid *inclusion_mask_grid = NULL,
  R3Point *points = NULL, int max_points = 0)
{
  // Parameters
  double max_world_distance = 2.0;
  double max_dd = max_world_distance * max_world_distance;

  // Initialize cluster
  cluster->parent = NULL;
  const RNScalar *position_x_grid_values = position_x_grid.GridValues();
  cluster->points.Insert(&position_x_grid_values[grid_index]);
  cluster->initial_segment_membership = 0;
  cluster->centroid = R3Point(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
  cluster->axes[0] = R3Vector(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
  cluster->axes[1] = R3Vector(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
  cluster->axes[2] = R3Vector(R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE, R2_GRID_UNKNOWN_VALUE);
  cluster->variances[0] = R2_GRID_UNKNOWN_VALUE;
  cluster->variances[1] = R2_GRID_UNKNOWN_VALUE;
  cluster->variances[2] = R2_GRID_UNKNOWN_VALUE;
  cluster->affinity = RN_EPSILON;

  // Check inclusion mask
  if (inclusion_mask_grid) {
    RNScalar inclusion_mask = inclusion_mask_grid->GridValue(grid_index);
    if ((inclusion_mask == 0) || (inclusion_mask == R2_GRID_UNKNOWN_VALUE)) return 1;
  }

  // Get initial segment membership
  int initial_segment_membership = 0;
  if (initial_segmentation_grid) {
    RNScalar initial_segment_membership_value = initial_segmentation_grid->GridValue(grid_index);
    if (initial_segment_membership_value == R2_GRID_UNKNOWN_VALUE) return 1;
    initial_segment_membership = (int) (initial_segment_membership_value + 0.5);
  }

  // Get point position
  RNScalar cx = position_x_grid.GridValue(grid_index);
  if (cx == R2_GRID_UNKNOWN_VALUE) return 1;
  RNScalar cy = position_y_grid.GridValue(grid_index);
  if (cy == R2_GRID_UNKNOWN_VALUE) return 1;
  RNScalar cz = position_z_grid.GridValue(grid_index);
  if (cz == R2_GRID_UNKNOWN_VALUE) return 1;
  R3Point position(cx, cy, cz);

  // Allocate temporary memory for neighbor points
  RNBoolean points_allocated = FALSE;
  if (!points) {
    int npoints = (2*pca_index_radius+1) * (2*pca_index_radius+1);
    points = new R3Point[npoints];
    points_allocated = TRUE;
  }

  // Find neighbor points
  int ci, cj, npoints = 0;
  position_x_grid.IndexToIndices(grid_index, ci, cj);
  for (int di = -pca_index_radius; di <= pca_index_radius; di++) {
    int i = ci + di;
    if ((i < 0) || (i >= position_x_grid.XResolution())) continue;
    for (int dj = -pca_index_radius; dj <= pca_index_radius; dj++) {
      int j = cj + dj;
      if ((j < 0) || (j >= position_x_grid.YResolution())) continue;
      if (inclusion_mask_grid) {
        if (inclusion_mask_grid->GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
      }
      if (initial_segmentation_grid && (initial_segment_membership > 0)) {
        RNScalar value = initial_segmentation_grid->GridValue(i, j);
        if (value == R2_GRID_UNKNOWN_VALUE) continue;
        if (!((int) (value + 0.5) & initial_segment_membership)) continue;
      }
      RNScalar x = position_x_grid.GridValue(i, j);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar y = position_y_grid.GridValue(i, j);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar z = position_z_grid.GridValue(i, j);
      if (z == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point p(x, y, z);
      RNLength dd = R3SquaredDistance(p, position);
      if (dd > max_dd) continue;
      if (npoints >= max_points) break;
      points[npoints++] = p;
    }
  }
      
  // Compute principle axes of neighbor points
  RNScalar variances[3];
  R3Point centroid = R3Centroid(npoints, points);
  R3Triad triad = R3PrincipleAxes(centroid, npoints, points, NULL, variances);

  // Delete points
  if (points && points_allocated) delete [] points;

  // Check variances to see if point should be included
  if ((pca_min_lambda1 > 0) && (variances[0] < pca_min_lambda1)) return 1; 
  if ((pca_min_lambda2 > 0) && (variances[1] < pca_min_lambda2)) return 1; 
  if ((pca_max_lambda2 > 0) && (variances[1] > pca_max_lambda2)) return 1; 
  if ((pca_max_lambda3 > 0) && (variances[2] > pca_max_lambda3)) return 1; 
  if ((pca_max_lambda21 > 0) && (variances[0] > 0) && (variances[1] / variances[0] > pca_max_lambda21)) return 1; 
  if ((pca_max_lambda31 > 0) && (variances[0] > 0) && (variances[2] / variances[0] > pca_max_lambda31)) return 1; 
  if ((pca_max_lambda32 > 0) && (variances[1] > 0) && (variances[2] / variances[1] > pca_max_lambda32)) return 1; 

  // POINT IS INCLUDED

  // Assign cluster centroid
  cluster->centroid = position;

  // Assign cluster axes
  cluster->axes[0] = triad[0];
  cluster->axes[1] = triad[1];
  cluster->axes[2] = triad[2];
 
  // Assign initial segment membership
  cluster->initial_segment_membership = initial_segment_membership;

  // Check seed mask
  if (seed_mask_grid) {
    RNScalar seed_mask = seed_mask_grid->GridValue(grid_index);
    if ((seed_mask == 0) || (seed_mask == R2_GRID_UNKNOWN_VALUE)) return 1;
  }

  // POINT IS SEED

  // Assign cluster variances
  cluster->variances[0] = variances[0];
  cluster->variances[1] = variances[1];
  cluster->variances[2] = variances[2];

  // Return success
  return 1;
}



static int
InitializeClusters(GSVScan *scan, Cluster *clusters, int nclusters,
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  const R2Grid& normal_x_grid, const R2Grid& normal_y_grid, const R2Grid& normal_z_grid,
  const R2Grid *initial_segmentation_grid = NULL, const R2Grid *seed_mask_grid = NULL, const R2Grid *inclusion_mask_grid = NULL)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int total_count = 0;
  int inclusion_candidate_count = 0;
  int seed_candidate_count = 0;
  int inclusion_count = 0;
  int seed_count = 0;

  // Parameters
  assert(nclusters == position_x_grid.NEntries());

  // Allocate temporary memory for points
  int npoints = (2*pca_index_radius+1) * (2*pca_index_radius+1);
  R3Point *points = new R3Point[npoints];

  // For every point
  for (int index = 0; index < nclusters; index++) {
    // Initialize cluster
    if (!InitializeCluster(scan, &clusters[index], index,
      position_x_grid, position_y_grid, position_z_grid,
      normal_x_grid, normal_y_grid, normal_z_grid,
      initial_segmentation_grid, seed_mask_grid, inclusion_mask_grid, 
      points, npoints)) {
      delete [] points;
      return 0;
    }

    // Gather some stats
    if (!inclusion_mask_grid || (inclusion_mask_grid->GridValue(index) > 0)) inclusion_candidate_count++;
    if (!seed_mask_grid || (seed_mask_grid->GridValue(index) > 0)) seed_candidate_count++;
    if (clusters[index].centroid[0] != R2_GRID_UNKNOWN_VALUE) inclusion_count++;
    if (clusters[index].variances[0] != R2_GRID_UNKNOWN_VALUE) seed_count++;
    total_count++;
  }

  // Delete points
  delete [] points;

  // Print debug message
  if (print_debug) printf("HEREA %d %d (%d) %d (%d) : %.3f s\n", total_count, 
    inclusion_count, inclusion_candidate_count, seed_count, seed_candidate_count, start_time.Elapsed());

  // Return success
  return 1;
}



static int
MergeClusters(GSVScan *scan, Cluster *clusters, int nclusters, 
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Create heap of cluster pairs 
  ClusterPair tmp;
  RNHeap<ClusterPair *> heap(&tmp, &tmp.affinity, &tmp.heapentry, FALSE);
  for (int index0 = 0; index0 < nclusters; index0++) {
    Cluster *cluster0 = &clusters[index0];
    if (cluster0->centroid[0] == R2_GRID_UNKNOWN_VALUE) continue;
    if (cluster0->centroid[1] == R2_GRID_UNKNOWN_VALUE) continue;
    if (cluster0->centroid[2] == R2_GRID_UNKNOWN_VALUE) continue;
    if (cluster0->variances[0] == R2_GRID_UNKNOWN_VALUE) continue;
    if (cluster0->variances[1] == R2_GRID_UNKNOWN_VALUE) continue;
    if (cluster0->variances[2] == R2_GRID_UNKNOWN_VALUE) continue;

    // Create pairs with neighbors
    int ix0, iy0;
    position_x_grid.IndexToIndices(index0, ix0, iy0);
    for (int s = -max_neighbor_scanline_index_difference; s <= max_neighbor_scanline_index_difference; s = (s<0) ? s/2 : ((s>0) ? s*2 : 1)) {
      for (int t = -max_neighbor_point_index_difference; t <= max_neighbor_point_index_difference; t++) {
        int index1;
        if ((s == 0) && (t == 0)) continue;
        int ix1 = ix0 + s;
        int iy1 = iy0 + t;
        if ((ix1 < 0) || (ix1 >= position_x_grid.XResolution())) continue;
        if ((iy1 < 0) || (iy1 >= position_x_grid.YResolution())) continue;
        position_x_grid.IndicesToIndex(ix1, iy1, index1);
        
        // Get cluster1
        Cluster *cluster1 = &clusters[index1];

        // Check cluster1
        if (cluster1->centroid[0] == R2_GRID_UNKNOWN_VALUE) continue;
        if (cluster1->centroid[1] == R2_GRID_UNKNOWN_VALUE) continue;
        if (cluster1->centroid[2] == R2_GRID_UNKNOWN_VALUE) continue;

        // Check if already considered this pair
        if (cluster1->variances[0] != R2_GRID_UNKNOWN_VALUE) {
          if (index0 >= index1) continue;
        }

        // Compute neighbor affinity
        RNScalar neighbor_affinity = NeighborAffinity(cluster0, cluster1);
        if (neighbor_affinity < 0) continue;
        if ((min_neighbor_affinity > 0) && (neighbor_affinity < min_neighbor_affinity)) continue;

        // Compute cluster affinity
        RNScalar pair_affinity = ClusterAffinity(cluster0, cluster1);
        if (pair_affinity < 0) continue;
        if ((min_pair_affinity > 0) && (pair_affinity < min_pair_affinity)) continue;

        // Create pair
        ClusterPair *pair = new ClusterPair();
        pair->clusters[0] = cluster0;
        pair->clusters[1] = cluster1;
        pair->affinity = pair_affinity;
        pair->heapentry = NULL;

        // Add pair to clusters
        cluster0->pairs.Insert(pair);
        cluster1->pairs.Insert(pair);

        // Insert pair into heap
        heap.Push(pair);
      }
    }
  }

  if (print_debug) {
    printf("HEREB %d (%.3f s)\n", heap.NEntries(), start_time.Elapsed());
    start_time.Read();
  }

  // Merge clusters hierarchically
  int pop_count = 0;
  int merge_count = 0;
  while (!heap.IsEmpty()) {
    // Get pair
    ClusterPair *pair01 = heap.Pop();
    pop_count++;

    // Get clusters
    Cluster *cluster0 = pair01->clusters[0];
    Cluster *cluster1 = pair01->clusters[1];

    // Check if clusters are OK 
    assert(cluster0->centroid[0] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster0->centroid[1] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster0->centroid[2] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster0->variances[0] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster0->variances[1] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster0->variances[2] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster1->centroid[0] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster1->centroid[1] != R2_GRID_UNKNOWN_VALUE);
    assert(cluster1->centroid[2] != R2_GRID_UNKNOWN_VALUE);

    // Check if affinity is high enough
    if (pair01->affinity > min_pair_affinity) {
      // Check if clusters have already been merged
      if (cluster0->parent || cluster1->parent) {
        // Find clusters
        Cluster *ancestor0 = cluster0;
        Cluster *ancestor1 = cluster1;
        while (ancestor0->parent) ancestor0 = ancestor0->parent;
        while (ancestor1->parent) ancestor1 = ancestor1->parent;
        if (ancestor0 != ancestor1) {
          // Find pair
          ClusterPair *pair = NULL;
          Cluster *ancestorA = (ancestor0->pairs.NEntries() < ancestor1->pairs.NEntries()) ? ancestor0 : ancestor1;
          Cluster *ancestorB = (ancestor0->pairs.NEntries() < ancestor1->pairs.NEntries()) ? ancestor1 : ancestor0;
          for (int j = 0; j < ancestorA->pairs.NEntries(); j++) {
            ClusterPair *tmp = ancestorA->pairs.Kth(j);
            if (tmp->clusters[0] == ancestorB) { pair = tmp; break; }
            if (tmp->clusters[1] == ancestorB) { pair = tmp; break; }
          }     
        
          // Create pair
          if (!pair) {
            RNScalar affinity = ClusterAffinity(ancestor0, ancestor1);
            if (affinity >= min_pair_affinity) {
              RNScalar denom = (ancestor0->affinity > ancestor1->affinity) ? ancestor0->affinity : ancestor1->affinity;
              RNScalar affinity_ratio = (denom > 0) ? affinity / denom : RN_INFINITY;
              if (affinity_ratio >= min_pair_affinity_ratio) {
                pair = new ClusterPair();
                pair->clusters[0] = ancestor0;
                pair->clusters[1] = ancestor1;
                pair->affinity = affinity;
                pair->heapentry = NULL;
                ancestor0->pairs.Insert(pair);
                ancestor1->pairs.Insert(pair);
                heap.Push(pair);
              }
            }
          }
        }
      }
      else {
        // Compute info for merged cluster
        int npoints0 = cluster0->points.NEntries();
        int npoints1 = cluster1->points.NEntries();
        RNScalar npoints = npoints0 + npoints1;
        RNScalar t0 = npoints0 / npoints;
        RNScalar t1 = npoints1 / npoints;
        R3Point centroid = t0 * cluster0->centroid + t1 * cluster1->centroid;
        // pair01->affinity = ClusterAffinity(cluster0, cluster1);

        // Print debug message
        if (print_debug) {
          if (0 && ((merge_count % 1000) == 0)) {
            RNScalar affinity = ClusterAffinity(cluster0, cluster1, print_debug);
            printf("  M %12.6f %12.6f : %15d %15d %15d : %9d %9d\n", 
                   pair01->affinity, affinity,
                   pop_count, merge_count, heap.NEntries(), 
                   npoints0, npoints1);
          }
        }
          
        // Merge clusters 
        Cluster *clusterA = (cluster0->points.NEntries() >= cluster1->points.NEntries()) ? cluster0 : cluster1;
        Cluster *clusterB = (cluster0->points.NEntries() >= cluster1->points.NEntries()) ? cluster1 : cluster0;
        clusterA->points.Append(clusterB->points);
        clusterA->centroid = centroid;
        clusterA->affinity = pair01->affinity;
        UpdateClusterTriad(clusterA, position_x_grid, position_y_grid, position_z_grid);
        assert(clusterA->variances[0] != R2_GRID_UNKNOWN_VALUE);
        clusterB->points.Empty(TRUE);
        clusterB->parent = clusterA;
        clusterB->affinity = -1;
          
        // Update statistics
        merge_count++;
      }
    }

    // Delete pair
    cluster0->pairs.Remove(pair01);
    cluster1->pairs.Remove(pair01);
    delete pair01;
  }

  // Print debug message
  if (print_debug) {
    printf("HEREC %d %d %d (%.3f s)\n", 
      pop_count, merge_count, heap.NEntries(),
      start_time.Elapsed());
  }

  // Return success
  return 1;
}



static int
ReassignClusters(GSVScan *scan, Cluster *clusters, int nclusters, RNArray<Cluster *>& sorted_clusters,
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  const R2Grid& normal_x_grid, const R2Grid& normal_y_grid, const R2Grid& normal_z_grid,
  const R2Grid *seed_mask_grid = NULL, const R2Grid *inclusion_mask_grid = NULL,
  const R2Grid *initial_segmentation_grid = NULL)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int nchanges = 0;

  // Make sorted list of merged clusters
  for (int i = 0; i < nclusters; i++) {
    Cluster *cluster = &clusters[i];
    if (cluster->parent) continue;
    sorted_clusters.Insert(cluster);
    // Update final cluster plane equation?
  }

  // Sort the list of clusters
  sorted_clusters.Sort(CompareClusters);

  // XXX TEMPORARY XX
  return 1;

  ////////////////////////////////////////////////////////////////////////

  // Create cluster index grid 
  R2Grid cluster_index_grid(position_x_grid.XResolution(), position_x_grid.YResolution());
  cluster_index_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  const RNScalar *position_x_grid_values = position_x_grid.GridValues();
  for (int i = 0; i < sorted_clusters.NEntries(); i++) {
    Cluster *cluster = sorted_clusters.Kth(i);
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      const RNScalar *point = cluster->points.Kth(j);
      int point_index = point - position_x_grid_values;
      cluster_index_grid.SetGridValue(point_index, i);
    }
  }

  // Allocate temporary memory
  int npoints = (2*pca_index_radius+1) * (2*pca_index_radius+1);
  R3Point *points = new R3Point[npoints];
  int *marks = new int [ sorted_clusters.NEntries() ];
  for (int i = 0; i < sorted_clusters.NEntries(); i++) marks[i] = 0;
  int current_mark = 1;

  // Make local adjustments to cluster index grid
  for (int cluster_index0 = sorted_clusters.NEntries()-1; cluster_index0 > 0; cluster_index0--) {
    Cluster *cluster0 = sorted_clusters.Kth(cluster_index0);
    for (int j = 0; j < cluster0->points.NEntries(); j++) {
      int ix0, iy0;
      const RNScalar *point0 = cluster0->points.Kth(j);
      int point_index0 = point0 - position_x_grid_values;
      if (position_x_grid.GridValue(point_index0) == R2_GRID_UNKNOWN_VALUE) continue;
      position_x_grid.IndexToIndices(point_index0, ix0, iy0);
      current_mark++;

      // Create cluster with just this one point
      Cluster standalone_cluster;
      if (!InitializeCluster(scan, &standalone_cluster, point_index0,
         position_x_grid, position_y_grid, position_z_grid,
         normal_x_grid, normal_y_grid, normal_z_grid,
         initial_segmentation_grid, NULL, inclusion_mask_grid, 
         points, npoints)) {
         continue;
      }

      // Check if there is a better nearby cluster
      // BEST_CLUSTER_INDEX DID NOT WORK AS LOCAL VARIABLE !!!???
      static int best_cluster_index; best_cluster_index = -1;  
      RNScalar best_affinity = min_neighbor_affinity;
      for (int s = -max_neighbor_scanline_index_difference; s <= max_neighbor_scanline_index_difference; s = (s<0) ? s/2 : ((s>0) ? s*2 : 1)) {
        for (int t = -max_neighbor_point_index_difference; t <= max_neighbor_point_index_difference; t++) {
          if ((s == 0) && (t == 0)) continue;
          int ix1 = ix0 + s;
          int iy1 = iy0 + t;
          if ((ix1 < 0) || (ix1 >= position_x_grid.XResolution())) continue;
          if ((iy1 < 0) || (iy1 >= position_x_grid.YResolution())) continue;
          RNScalar cluster_index1_value = cluster_index_grid.GridValue(ix1, iy1);
          if (cluster_index1_value == R2_GRID_UNKNOWN_VALUE) continue;
          int cluster_index1 = (int) (cluster_index1_value + 0.5);
          if (cluster_index1 >= cluster_index0) continue;
          if (marks[cluster_index1] == current_mark) continue;
          marks[cluster_index1] = current_mark;
          Cluster *cluster1 = sorted_clusters.Kth(cluster_index1);
          RNScalar affinity = ClusterAffinity(cluster1, &standalone_cluster);
          if (affinity > best_affinity) {
            best_cluster_index = cluster_index1;
            best_affinity = affinity;
          }
        }
      }

      // Check if found better cluster
      if ((best_cluster_index >= 0) && (best_cluster_index != cluster_index0)) {
        // Reassign point0 to best cluster
        cluster_index_grid.SetGridValue(point_index0, best_cluster_index);
        Cluster *best_cluster = sorted_clusters.Kth(best_cluster_index);
        best_cluster->points.Insert(point0);
        RNArrayEntry *entry = cluster0->points.KthEntry(j);
        cluster0->points.EntryContents(entry) = cluster0->points.Tail();
        cluster0->points.RemoveTail(); j--;
        nchanges++;
      }          
    }
  }

  // Delete temporary memory
  delete [] points;
  delete [] marks;

  // Print debug message
  if (print_debug) {
    printf("HERED %d (%.3f s)\n", 
      nchanges, start_time.Elapsed());
  }

  // Return success
  return 1;
}



static int
OrientClusters(GSVScan *scan, const RNArray<Cluster *>& clusters,
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  const R2Grid& normal_x_grid, const R2Grid& normal_y_grid, const R2Grid& normal_z_grid)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int nchanges = 0;

  // Orient plane for all clusters
  const RNScalar *position_x_grid_values = position_x_grid.GridValues();
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster = clusters.Kth(i);

    // Determine orientation 
    RNScalar sum = 0;
    R3Vector cluster_normal = cluster->axes[2];
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      const RNScalar *point = cluster->points.Kth(j);
      int grid_index = point - position_x_grid_values;
      RNCoord tx = normal_x_grid.GridValue(grid_index);
      if (tx == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord ty = normal_y_grid.GridValue(grid_index);
      if (ty == R2_GRID_UNKNOWN_VALUE) continue;
      RNCoord tz = normal_z_grid.GridValue(grid_index);
      if (tz == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector point_normal(tx, ty, tz);
      RNScalar dot = cluster_normal.Dot(point_normal);
      sum += dot;
    }

    // Flip if cluster normal is backfacing
    if (sum < 0) {
      cluster->axes[0].Flip();
      cluster->axes[1].Flip();
      cluster->axes[2].Flip();
      nchanges++;
    }
  }

  // Print debug message
  if (print_debug) {
    printf("HEREE %d (%.3f s)\n", 
      nchanges, start_time.Elapsed());
  }

  // Return success
  return 1;
}



static int
WritePlaneGrids(GSVScan *scan, const RNArray<Cluster *>& clusters,
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid,
  const char *output_image_directory, const char *output_image_basename, const char *image_parameterization,
  const R2Grid *seed_mask_grid = NULL, const R2Grid *inclusion_mask_grid = NULL)
{
  // Get convenient variables
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  char image_name[4096];

  // Create output grids
  R2Grid segment_id_grid(position_x_grid);     
  segment_id_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid segment_size_grid(segment_id_grid);
  R2Grid segment_A_grid(segment_id_grid);
  R2Grid segment_B_grid(segment_id_grid);
  R2Grid segment_C_grid(segment_id_grid);
  R2Grid segment_D_grid(segment_id_grid);

  // Fill output grids
  int cluster_count = 0;
  const RNScalar *position_x_grid_values = position_x_grid.GridValues();
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster = clusters.Kth(i);

    // Create plane equation
    R3Plane plane(cluster->centroid, cluster->axes[2]);

    if (print_debug) {
      if (cluster->points.NEntries() > 100) {
        printf("  C %6d : %6d ( %7.3f %7.3f %7.3f ) ( %7.3f %7.3f %7.3f %7.3f ) %7.3f\n", i, cluster->points.NEntries(), 
               cluster->centroid[0], cluster->centroid[1], cluster->centroid[2], 
               plane[0], plane[1], plane[2], plane[3], 
               cluster->affinity);
      }
    }

    // Update grids
    int point_count = 0;
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      const RNScalar *point = cluster->points.Kth(j);
      int point_index = point - position_x_grid_values;
      if (position_x_grid.GridValue(point_index) == R2_GRID_UNKNOWN_VALUE) {
        continue;
      }
      if (inclusion_mask_grid) {
        if (inclusion_mask_grid->GridValue(point_index) == 0) continue;
        if (inclusion_mask_grid->GridValue(point_index) == R2_GRID_UNKNOWN_VALUE) continue;
      }
      segment_id_grid.SetGridValue(point_index, cluster_count);
      segment_size_grid.SetGridValue(point_index, cluster->points.NEntries());
      segment_A_grid.SetGridValue(point_index, plane.A());
      segment_B_grid.SetGridValue(point_index, plane.B());
      segment_C_grid.SetGridValue(point_index, plane.C());
      segment_D_grid.SetGridValue(point_index, plane.D());
      point_count++;
    }

    // Update cluster counter
    if (point_count > 0) cluster_count++;
  }


  ////////////////////////////////////////////////////////////////////////

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sId.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_id_grid.Write(image_name);

  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sSize.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_size_grid.Write(image_name);

  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sA.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_A_grid.Write(image_name);                                                                                                                        
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sB.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_B_grid.Write(image_name);                                                                                                                        
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sC.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_C_grid.Write(image_name);                                                                                                                        
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sD.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  segment_D_grid.Write(image_name);

  // Return success
  return 1;
}



static int
SegmentPlanes(GSVScan *scan, 
  const char *input_image_directory, const char *input_initial_segmentation_basename,  
  const char *input_seed_mask_basename, const char *input_inclusion_mask_basename,
  const char *output_image_directory, const char *output_image_basename, 
  const char *image_parameterization)
{
  // Get convenient variables
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  char image_name[4096];

  ////////////////////////////////////////////////////////////////////////

  // Check if already computed
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sId.grd", output_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, output_image_basename);
  if (RNFileExists(image_name)) return 1;

  ////////////////////////////////////////////////////////////////////////

  // Read grids
  R2Grid position_x_grid;
  R2Grid position_y_grid;
  R2Grid position_z_grid;
  R2Grid normal_x_grid;
  R2Grid normal_y_grid;
  R2Grid normal_z_grid;
  R2Grid viewpoint_depth_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_z_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normal_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalY.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normal_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalZ.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normal_z_grid.Read(image_name)) return 0;

  // Read initial segmentation grid
  R2Grid *initial_segmentation_grid = NULL;
  if (input_initial_segmentation_basename) {
    sprintf(image_name, "%s/%s/%02d_%02d_%s_%s.grd", input_image_directory, run->Name(), segment_index, scan_index, 
      image_parameterization, input_initial_segmentation_basename);
    initial_segmentation_grid = new R2Grid();
    if (!initial_segmentation_grid->Read(image_name)) return 0;
  }

  // Create inclusion mask
  R2Grid inclusion_mask(position_x_grid);
  inclusion_mask.Threshold(-FLT_MAX, R2_GRID_UNKNOWN_VALUE, 1);
  if (mask_sensor_vehicle) {
    MaskPointsOnSensorVehicle(inclusion_mask);
  }
  if (mask_shadow_boundaries_from_inclusion) {
    R2Grid boundary_type_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
    if (!boundary_type_grid.Read(image_name)) return 0;
    MaskBoundaries(inclusion_mask, boundary_type_grid, GSV_SHADOW_BOUNDARIES);
  }
  if (input_inclusion_mask_basename) {
    R2Grid mask_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_%s_%s.grd", input_image_directory, run->Name(), segment_index, scan_index, 
      image_parameterization, input_inclusion_mask_basename);
    if (!mask_grid.Read(image_name)) return 0;
    mask_grid.Substitute(0, R2_GRID_UNKNOWN_VALUE);
    mask_grid.Threshold(RN_EPSILON, R2_GRID_UNKNOWN_VALUE, R2_GRID_KEEP_VALUE); // This is a hack for convenience of curbs, etc.
    inclusion_mask.Mask(mask_grid);
  }
  
  // Create seed mask
  R2Grid seed_mask(inclusion_mask);
  if (mask_inconsistent_normals_from_seeds) {
    MaskInconsistentNormals(seed_mask, normal_x_grid, normal_y_grid, normal_z_grid, max_pair_axis_angle[2]); 
  }
  if (mask_boundaries_from_seeds) {
    R2Grid boundary_type_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
    if (!boundary_type_grid.Read(image_name)) return 0;
    MaskBoundaries(seed_mask, boundary_type_grid, GSV_SHADOW_BOUNDARIES | GSV_SILHOUETTE_BOUNDARIES);
    // R2Grid viewpoint_depth_grid;
    // sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointDepth.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
    // if (!viewpoint_depth_grid.Read(image_name)) return 0;
    // MaskOcclusionBoundaries(seed_mask, viewpoint_depth_grid);
  }
  if (input_seed_mask_basename) {
    R2Grid mask_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_%s_%s.grd", input_image_directory, run->Name(), segment_index, scan_index, 
      image_parameterization, input_seed_mask_basename);
    if (!mask_grid.Read(image_name)) return 0;
    mask_grid.Substitute(0, R2_GRID_UNKNOWN_VALUE);
    mask_grid.Threshold(RN_EPSILON, R2_GRID_UNKNOWN_VALUE, R2_GRID_KEEP_VALUE); // This is a hack for convenience of curbs, etc.
    seed_mask.Mask(mask_grid);
  }

  // Allocate clusters
  int nclusters = position_x_grid.NEntries();
  Cluster *clusters = new Cluster [ nclusters ];
  if (!clusters) {
    fprintf(stderr, "Unable to allocate clusters\n");
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    return 0;
  }

#if 0
  // Write debug grids
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%s.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization, "InclusionMask");
  inclusion_mask.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%s.grd", input_image_directory, run->Name(), segment_index, scan_index, image_parameterization, "SeedMask");
  seed_mask.Write(image_name);
#endif
  
  // Initialize clusters
  if (!InitializeClusters(scan, clusters, nclusters,
    position_x_grid, position_y_grid, position_z_grid,
    normal_x_grid, normal_y_grid, normal_z_grid,
    initial_segmentation_grid, &seed_mask, &inclusion_mask)) {
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    delete [] clusters;
    return 0;
  }

  // Merge clusters
  if (!MergeClusters(scan, clusters, nclusters,
    position_x_grid, position_y_grid, position_z_grid)) {
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    delete [] clusters;
    return 0;
  }

  // Reassign points to clusters
  RNArray<Cluster *> sorted_clusters;
  if (!ReassignClusters(scan, clusters, nclusters, sorted_clusters,
    position_x_grid, position_y_grid, position_z_grid,
    normal_x_grid, normal_y_grid, normal_z_grid,
    initial_segmentation_grid, &seed_mask, &inclusion_mask)) {
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    delete [] clusters;
    return 0;
  }

  // Orient clusters
  if (!OrientClusters(scan, sorted_clusters,
    position_x_grid, position_y_grid, position_z_grid,
    normal_x_grid, normal_y_grid, normal_z_grid)) {
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    delete [] clusters;
    return 0;
  }

  // Write clusters to grids
  if (!WritePlaneGrids(scan, sorted_clusters,
    position_x_grid, position_y_grid, position_z_grid,
    output_image_directory, output_image_basename, image_parameterization, 
    &seed_mask, &inclusion_mask)) {
    if (initial_segmentation_grid) delete initial_segmentation_grid;
    delete [] clusters;
    return 0;
  }

  // Delete initial segmentation grid
  if (initial_segmentation_grid) delete initial_segmentation_grid;

  // Delete clusters
  delete [] clusters;

  // Print statistics
  if (print_verbose) {
    printf("  Segmented %s images for %s %02d %02d\n", image_parameterization, run->Name(), segment_index, scan_index);
    fflush(stdout);
  }

  // Return success
  return 1;
}


static int
SegmentPlanes(GSVScene *scene, 
  const char *input_image_directory, const char *input_initial_segmentation_basename, 
  const char *input_seed_mask_basename, const char *input_inclusion_mask_basename, 
  const char *output_image_directory, const char *output_image_basename, 
  const char *image_parameterization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Segmenting %s images ...\n", image_parameterization);
    fflush(stdout);
  }

  // Segment images for all scans
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (!SegmentPlanes(scan, input_image_directory, input_initial_segmentation_basename, 
          input_seed_mask_basename, input_inclusion_mask_basename,
          output_image_directory, output_image_basename, image_parameterization)) return 0;
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

static int 
ParseArgs(int argc, char **argv)
{
  // By default ...
  int segment_default_images = 1;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-debug")) print_debug = 1; 
      else if (!strcmp(*argv, "-planar_grids")) write_planar_grids = 1; 
      else if (!strcmp(*argv, "-xy_spans")) write_xy_spans = 1; 
      else if (!strcmp(*argv, "-seed_mask_name")) { argc--; argv++; input_seed_mask_basename = *argv; }
      else if (!strcmp(*argv, "-inclusion_mask_name")) { argc--; argv++; input_inclusion_mask_basename = *argv; }
      else if (!strcmp(*argv, "-initial_segmentation_name")) { argc--; argv++; input_initial_segmentation_basename = *argv; }
      else if (!strcmp(*argv, "-output_name")) { argc--; argv++; output_image_basename = *argv; }
      else if (!strcmp(*argv, "-SA")) { segment_SA_images = 1; segment_default_images = 0; }
      else if (!strcmp(*argv, "-DA")) { segment_DA_images = 1; segment_default_images = 0; }
      else if (!strcmp(*argv, "-DH")) { segment_DH_images = 1; segment_default_images = 0; }
      else if (!strcmp(*argv, "-small_planes")) { 
        output_image_basename = "SmallPlane";
        max_neighbor_scanline_index_difference = 32;
        max_neighbor_point_index_difference = 8;
        max_neighbor_centroid_distance = 2;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 0.1;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = RN_PI / 8.0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 0.1;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = RN_PI / 64.0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = RN_PI / 32.0;
      }
      else if (!strcmp(*argv, "-very_tight_planes")) { 
        output_image_basename = "VeryTightPlane";
        max_neighbor_scanline_index_difference = 32;
        max_neighbor_point_index_difference = 8;
        max_neighbor_centroid_distance = 2;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 0.25;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = RN_PI / 8.0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 0.25;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = RN_PI / 32.0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = RN_PI / 32.0;
      }
      else if (!strcmp(*argv, "-tight_planes")) { 
        output_image_basename = "TightPlane";
        max_neighbor_scanline_index_difference = 32;
        max_neighbor_point_index_difference = 8;
        max_neighbor_centroid_distance = 2;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 0.5;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = RN_PI / 8.0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 0.5;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = RN_PI / 32.0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = RN_PI / 16.0;
      }
      else if (!strcmp(*argv, "-large_planes")) { 
        output_image_basename = "LargePlane";
        max_neighbor_scanline_index_difference = 64;
        max_neighbor_point_index_difference = 8;
        max_neighbor_centroid_distance = 4;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 2;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = RN_PI / 4.0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 1;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = RN_PI / 8.0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = RN_PI / 8.0;
      }
      else if (!strcmp(*argv, "-loose_planes")) { 
        output_image_basename = "LoosePlane";
        max_neighbor_scanline_index_difference = 64;
        max_neighbor_point_index_difference = 8;
        max_neighbor_centroid_distance = 4;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 2;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = RN_PI / 3.0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 2;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = RN_PI / 8.0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = RN_PI / 8.0;
      }
      else if (!strcmp(*argv, "-boundary_lines")) { 
        input_initial_segmentation_basename = "BoundaryType";
        input_inclusion_mask_basename = "BoundaryType";
        output_image_basename = "BoundaryLine";
        pca_max_lambda2 = 0.01;
        pca_max_lambda21 = 0.25;
        mask_shadow_boundaries_from_inclusion = 1;
        mask_boundaries_from_seeds = 0;
        mask_inconsistent_normals_from_seeds = 0;
        max_neighbor_scanline_index_difference = 4;
        max_neighbor_point_index_difference = 4;
        max_neighbor_centroid_distance = 4;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0.25;
        max_neighbor_offplane_distance[2] = 0.25;
        max_neighbor_axis_angle[0] = RN_PI / 8.0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = 0;
        max_pair_centroid_distance = 100;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0.25;
        max_pair_offplane_distance[2] = 0.25;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = RN_PI / 8.0;
        max_pair_offplane_angle[2] = RN_PI / 8.0;
        max_pair_axis_angle[0] = RN_PI / 8.0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = 0;
      }
      else if (!strcmp(*argv, "-small_objects")) { 
        output_image_basename = "SmallObject";
        mask_boundaries_from_seeds = 0;
        mask_inconsistent_normals_from_seeds = 0;
        pca_min_lambda1 = 0;
        pca_min_lambda2 = 0;
        pca_max_lambda2 = 0;
        pca_max_lambda3 = 0;
        pca_max_lambda21 = 0;
        pca_max_lambda31 = 0;
        pca_max_lambda32 = 0;
        max_neighbor_scanline_index_difference = 8;
        max_neighbor_point_index_difference = 4;
        max_neighbor_centroid_distance = 0.5;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 0;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = 0;
        max_pair_centroid_distance = 8;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 0;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = 0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = 0;
        min_pair_affinity_ratio = 0.1;
      }
      else if (!strcmp(*argv, "-superpixels")) { 
        output_image_basename = "Superpixels";
        mask_boundaries_from_seeds = 0;
        mask_inconsistent_normals_from_seeds = 0;
        pca_min_lambda1 = 0;
        pca_min_lambda2 = 0;
        pca_max_lambda2 = 0;
        pca_max_lambda3 = 0;
        pca_max_lambda21 = 0;
        pca_max_lambda31 = 0;
        pca_max_lambda32 = 0;
        max_neighbor_scanline_index_difference = 8;
        max_neighbor_point_index_difference = 4;
        max_neighbor_centroid_distance = 0.5;
        max_neighbor_offplane_distance[0] = 0;
        max_neighbor_offplane_distance[1] = 0;
        max_neighbor_offplane_distance[2] = 0;
        max_neighbor_axis_angle[0] = 0;
        max_neighbor_axis_angle[1] = 0;
        max_neighbor_axis_angle[2] = 0;
        max_pair_centroid_distance = 8;
        max_pair_offplane_distance[0] = 0;
        max_pair_offplane_distance[1] = 0;
        max_pair_offplane_distance[2] = 0;
        max_pair_offplane_angle[0] = 0;
        max_pair_offplane_angle[1] = 0;
        max_pair_offplane_angle[2] = 0;
        max_pair_axis_angle[0] = 0;
        max_pair_axis_angle[1] = 0;
        max_pair_axis_angle[2] = 0;
        min_pair_affinity_ratio = 0.5;
      }
      else if (!strcmp(*argv, "-min_pair_affinity_ratio")) { 
        argc--; argv++; min_pair_affinity_ratio = atof(*argv); 
      }
      else if (!strcmp(*argv, "-min_pair_affinity")) { 
        argc--; argv++; min_pair_affinity = atof(*argv); 
      }
      else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!input_image_directory) input_image_directory = *argv;
      else if (!output_image_directory) output_image_directory = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsvsegment input_scene [input_image_directory] [output_image_directory] [options]\n");
    return FALSE;
  }

  // Check input image directory
  if (!input_image_directory) {
    input_image_directory = "gsv_data/laser_images";
  }

  // Check output image directory
  if (!output_image_directory) {
    output_image_directory = "gsv_data/laser_images";
  }

  // Check image selections
  if (segment_default_images) segment_DA_images = 1;

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
    
  // Segment SA images
  if (segment_SA_images) {
    if (!SegmentPlanes(scene, input_image_directory, input_initial_segmentation_basename, 
      input_seed_mask_basename, input_inclusion_mask_basename,
      output_image_directory, output_image_basename, "SA")) exit(-1);
    if (write_planar_grids) {
      if (!WritePlanarGrids(scene, input_image_directory, 
        output_image_directory, output_image_basename, "SA")) exit(-1);
    }
    if (write_xy_spans) {
      if (!WriteXYSpans(scene, input_image_directory, 
        output_image_directory, output_image_basename, "SA")) exit(-1);
    }
  }

  // Segment DA images
  if (segment_DA_images) {
    if (!SegmentPlanes(scene, input_image_directory, input_initial_segmentation_basename, 
      input_seed_mask_basename, input_inclusion_mask_basename,
      output_image_directory, output_image_basename, "DA")) exit(-1);
    if (write_planar_grids && !segment_SA_images) {
      if (!WritePlanarGrids(scene, input_image_directory, 
        output_image_directory, output_image_basename, "DA")) exit(-1);
    }
    if (write_xy_spans && !segment_SA_images) {
      if (!WriteXYSpans(scene, input_image_directory, 
        output_image_directory, output_image_basename, "DA")) exit(-1);
    }
  }

  // Segment DH images
  if (segment_DH_images) {
    if (!SegmentPlanes(scene, input_image_directory, input_initial_segmentation_basename, 
      input_seed_mask_basename, input_inclusion_mask_basename, 
      output_image_directory, output_image_basename, "DH")) exit(-1);
    if (write_planar_grids && !segment_SA_images && !segment_DA_images) {
      if (!WritePlanarGrids(scene, input_image_directory, 
        output_image_directory, output_image_basename, "DH")) exit(-1);
    }
    if (write_xy_spans && !segment_SA_images && !segment_DA_images) {
      if (!WriteXYSpans(scene, input_image_directory, 
        output_image_directory, output_image_basename, "DH")) exit(-1);
    }
  }

  // Return success 
  return 0;
}

















