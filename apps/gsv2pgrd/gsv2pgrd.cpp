// Source file for the google viewer program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"
#include "GSV/GSV.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static const char *input_scene_name = NULL;
static const char *input_image_directory = NULL;
static const char *output_image_directory = NULL;
static const char *image_basename = "Large";
static const char *image_parameterization = "SA";
static int horizontal = FALSE;
static int vertical = FALSE;
static int oblique = FALSE;
static int print_verbose = 0;
static int print_debug = 0;



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
  const char *input_image_directory, const char *output_image_directory, 
  const char *image_basename, const char *image_parameterization, 
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
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneA.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_a_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneB.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_b_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneC.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_c_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneD.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_d_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneId.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_id_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneSize.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_size_grid.Read(image_name)) return 0;

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

    // Check plane orientation
    RNScalar normal_tolerance = 0.01;
    if (RNIsEqual(normal.Z(), 0.0, normal_tolerance)) { if (!vertical) continue; }
    else if (RNIsEqual(normal.Z(), 1.0, normal_tolerance)) { if (!horizontal) continue; }
    else if (RNIsEqual(normal.Z(), -1.0, normal_tolerance)) { if (!horizontal) continue; }
    else { if (!oblique) continue; }

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
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p10, p11, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p10, p11, d00, d10, d11);
          height_planar_grid.RasterizeWorldTriangle(p00, p10, p11, h00, h10, h11);
          density_planar_grid.RasterizeWorldTriangle(p00, p10, p11, s00, s10, s11);
          overlap_planar_grid.RasterizeWorldTriangle(p00, p10, p11, 1.0);
        }
        if (p01[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d00 = R3SignedDistance(plane, p00);
          RNScalar d01 = R3SignedDistance(plane, p01);
          RNScalar d11 = R3SignedDistance(plane, p11);
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p11, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p11, p01, d00, d11, d01);
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
          RNScalar h00 = height_grid.GridValue(ix, iy);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar s00 = density_grid.GridValue(ix, iy);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p00, p10, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p00, p10, p01, d00, d10, d01);
          height_planar_grid.RasterizeWorldTriangle(p00, p10, p01, h00, h10, h01);
          density_planar_grid.RasterizeWorldTriangle(p00, p10, p01, s00, s10, s01);
          overlap_planar_grid.RasterizeWorldTriangle(p00, p10, p01, 1.0);
        }
        if (p11[0] != R2_GRID_UNKNOWN_VALUE) {
          RNScalar d01 = R3SignedDistance(plane, p01);
          RNScalar d10 = R3SignedDistance(plane, p10);
          RNScalar d11 = R3SignedDistance(plane, p11);
          RNScalar h10 = height_grid.GridValue(ix+1, iy);
          RNScalar h01 = height_grid.GridValue(ix, iy+1);
          RNScalar h11 = height_grid.GridValue(ix+1, iy+1);
          RNScalar s10 = density_grid.GridValue(ix+1, iy);
          RNScalar s01 = density_grid.GridValue(ix, iy+1);
          RNScalar s11 = density_grid.GridValue(ix+1, iy+1);
          index_planar_grid.RasterizeWorldTriangle(p10, p11, p01, grid_index);
          displacement_planar_grid.RasterizeWorldTriangle(p10, p11, p01, d10, d11, d01);
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
    height_planar_grid.Mask(mask);
    density_planar_grid.Mask(mask);
    index_planar_grid.Mask(mask);

    // Compute means
    displacement_planar_grid.Divide(overlap_planar_grid);
    height_planar_grid.Divide(overlap_planar_grid);
    density_planar_grid.Divide(overlap_planar_grid);
    index_planar_grid.Divide(overlap_planar_grid); // NOT RIGHT

    // Write planar grids
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Displacement.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!displacement_planar_grid.WriteFile(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_%sGridIndex.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id, image_parameterization);
    if (!index_planar_grid.WriteFile(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Height.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!height_planar_grid.WriteFile(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Density.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!density_planar_grid.WriteFile(image_name)) return 0;

    // Print debug message
    if (print_debug) {
      printf("  G %20s %2d %2d : %6.3f %6.3f %6.3f %6.3f : %6.3f %6.3f %6.3f : %d\n", 
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
WriteColorGrids(GSVScan *scan, 
  const char *input_image_directory, const char *output_image_directory, 
  const char *image_basename, const char *image_parameterization,
  RNScalar color_grid_spacing = 0.05,
  RNLength max_distance = 1000, RNAngle max_angle = RN_PI_OVER_TWO)
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
  RNScalar min_distance = 1.0;

  // Read plane segment id grid
  R2Grid plane_segment_id_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_%sPlaneId.grd", 
    input_image_directory, run->Name(), segment_index, scan_index, 
    image_parameterization, image_basename);
  if (!plane_segment_id_grid.Read(image_name)) return 0;
  int max_segment_id = (int) (plane_segment_id_grid.Maximum() + 0.5);

  // Create mesh
  GSVMesh *mesh = scan->Mesh();
  // if (!mesh) continue;
    
  // Process segments
  for (int id = 0; id <= max_segment_id; id++) {
    // Check if segment files exist
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Displacement.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!RNFileExists(image_name)) continue;

    // Read displacement grid
    R3PlanarGrid displacement_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Displacement.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, id);
    if (!displacement_grid.ReadFile(image_name)) continue;
    displacement_grid.Substitute(R2_GRID_UNKNOWN_VALUE, 0.0);

    // Consider tapestries relevant to scan
    int ntapestries[3] = { 6, 3, 6 };
    int tapestries[3][6] = { { 4, 5, 6, 7, 0, 8 }, { 0 }, { 0, 1, 2, 3, 4, 8 } };
    for (int k = 0; k < ntapestries[scan_index]; k++) {
      int it = tapestries[scan_index][k];
      GSVTapestry *tapestry = segment->Tapestry(it);
      GSVCamera *camera = tapestry->Camera();
      if (!camera) continue;
      
      // Consider each image of tapestry
      for (int ii = 0; ii < tapestry->NImages(); ii++) {
        GSVImage *image = tapestry->Image(ii);
        int image_index = image->PanoramaIndex();
        GSVPanorama *panorama = image->Panorama();
        int panorama_index = panorama->SegmentIndex();

        // Get image pose info
        const GSVPose& pose = image->Pose();
        const R3Point& viewpoint = pose.Viewpoint();
        const R3Vector& towards = pose.Towards();
        const R3Vector& up = pose.Up();
        RNScalar xfov = 0.5 * camera->XFov();
        RNScalar yfov = 0.5 * camera->YFov();
        R3Frustum frustum(viewpoint, towards, up, xfov, yfov, min_distance, max_distance);
        
        // Check if within viewing angle
        const R3Vector& normal = displacement_grid.Plane().Normal();
        RNAngle angle = R3InteriorAngle(-normal, towards);
        if (angle > RN_PI_OVER_TWO) continue;
        // if (angle > max_angle) continue;
        
        // Check if within distance range
        const R3Box& world_box = displacement_grid.WorldBox();
        RNLength distance = R3Distance(viewpoint, world_box);
        if ((max_distance > 0) && (distance > max_distance)) continue;

        // Check if within view frustum
        if (!frustum.Intersects(world_box)) continue;
        
        // Check if already done
        sprintf(image_name, "%s/%s/%02d_%02d_%06d_%06d_%02d_ST_Color.bmp", 
          output_image_directory, run->Name(), segment_index, scan_index, id, panorama_index, image_index);
        if (RNFileExists(image_name)) continue;

        // Get distorted image
        R2Image *distorted_image = image->DistortedImage();
        if (!distorted_image) continue;

        // Create distance grid
        RNScalar scale = 1;
        int xres = displacement_grid.XResolution();
        int yres = displacement_grid.YResolution();
        R3PlanarGrid distance_grid(displacement_grid); 
        distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
        // if (color_grid_spacing != RN_UNKNOWN) {
        //  scale = displacement_grid.WorldToGridScaleFactor() / color_grid_spacing;
        //  int xres = scale * displacement_grid.XResolution();
        //  int yres = scale * displacement_grid.YResolution();
        //  distance_grid.Resample(xres, yres);
        //}
        
        // Create RGB image
        R2Image rgb_image(xres, yres);
        RNScalar sum = 0.0;

        // Compute color for every pixel 
        for (int ix = 0; ix < distance_grid.XResolution(); ix++) {
          for (int iy = 0; iy < distance_grid.YResolution(); iy++) {
            // Get world position
            R3Point world_position = distance_grid.WorldPosition(ix, iy);
            RNScalar displacement = displacement_grid.GridValue(ix / scale, iy / scale);
            if (displacement == R2_GRID_UNKNOWN_VALUE) continue;
            world_position += displacement * distance_grid.Plane().Normal();

            // Check if within distance range
            RNLength distance = R3Distance(viewpoint, world_position);
            if ((max_distance > 0) && (distance > max_distance)) continue;
            
            // Check if closest so far
            RNScalar old_distance = distance_grid.GridValue(ix, iy);
            if ((old_distance != R2_GRID_UNKNOWN_VALUE) && (distance > old_distance)) continue;

            // Check if within view frustum
            if (!frustum.Intersects(world_position)) continue;
#if 1
            // Check if closest intersection
            R3Ray ray(viewpoint, world_position);
            R3MeshIntersection intersection;
            const RNLength visibility_tolerance = 1.0;
            if (mesh && mesh->Intersection(ray, &intersection, 0, distance)) {
              if (intersection.t < distance - visibility_tolerance) continue;
            }
#endif
            // Retrieve color of point from distorted image
            R2Point image_position = image->DistortedPosition(world_position);
            if (image_position.X() < 0) continue;
            if (image_position.Y() < 0) continue;
            int image_ix1 = (int) image_position.X();
            int image_iy1 = (int) image_position.Y();
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
            RNRgb colorA = (1-tx)*color11 + tx*color21;
            RNRgb colorB = (1-tx)*color12 + tx*color22;
            RNRgb color = (1-ty)*colorA + ty*colorB; 

            // Update images
            rgb_image.SetPixelRGB(ix, iy, color);
            distance_grid.SetGridValue(ix, iy, distance);
            sum += color.Luminance();
          }
        }

        // Write color image
        if (sum > 0.0) {
          sprintf(image_name, "%s/%s/%02d_%02d_%06d_%06d_%02d_ST_Color.bmp", 
            output_image_directory, run->Name(), segment_index, scan_index, id, panorama_index, image_index);
          rgb_image.Write(image_name);

          // Print debug message
          if (print_debug) {
            printf("  C %20s %2d %2d %6d %6d %2d\n", run->Name(), segment_index, scan_index, id, panorama_index, image_index);
          }
        }

        // Delete distorted image
        delete distorted_image;
      }
    }
  }

  // Delete mesh
  if (mesh) delete mesh;

  // Return success
  return 1;
}




static int
WritePlanarGrids(GSVScene *scene, 
  const char *input_image_directory, 
  const char *output_image_directory, const char *image_basename, 
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
          output_image_directory, image_basename, 
          image_parameterization)) return 0;
        if (!WriteColorGrids(scan, input_image_directory, 
          output_image_directory, image_basename, 
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
      else if (!strcmp(*argv, "-debug")) print_debug = 1; 
      else if (!strcmp(*argv, "-horizontal")) horizontal = TRUE; 
      else if (!strcmp(*argv, "-vertical")) vertical = TRUE; 
      else if (!strcmp(*argv, "-oblique"))  oblique = TRUE; 
      else if (!strcmp(*argv, "-name")) { argc--; argv++; image_basename = *argv; }
      else if (!strcmp(*argv, "-parameterization")) { argc--; argv++; image_parameterization = *argv; }
      else if (!strcmp(*argv, "-SA")) { image_parameterization = "SA"; }
      else if (!strcmp(*argv, "-DA")) { image_parameterization = "DA"; }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
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
    fprintf(stderr, "Usage: gsv2pgrd input_scene [input_image_directory] [output_image_directory] [options]\n");
    return FALSE;
  }

  // Check input image directory
  if (!input_image_directory) {
    input_image_directory = "gsv_data/laser_images";
  }

  // Check output image directory
  if (!output_image_directory) {
    output_image_directory = "gsv_data/planar_images";
  }

  // Check plane types
  if (!horizontal && !vertical && !oblique) {
    horizontal = vertical = oblique = TRUE;
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
    
  // Create output images
  if (!WritePlanarGrids(scene, 
    input_image_directory, output_image_directory, 
    image_basename, image_parameterization)) exit(-1);

  // Return success 
  return 0;
}

















