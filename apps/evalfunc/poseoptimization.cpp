// Source file for GSV pose optimization



////////////////////////////////////////////////////////////////////////
// Select solvers (guides compilation in RNMath/RNPolynomial.h)
////////////////////////////////////////////////////////////////////////

/*#define RN_USE_CSPARSE
#define RN_USE_SPLM
#define RN_USE_CERES
#define RN_USE_MINPACK
*/


////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "RNMath/RNMath.h"
#include "R3CatmullRomSpline.h"
#include "poseoptimization.h"



/////////////////////////////////////////////////////////////////////////
// GSVDescriptor utility functions
////////////////////////////////////////////////////////////////////////

GSVDescriptor::
GSVDescriptor(int descriptor_type, RNScalar *values, int nvalues)
  : descriptor_type(descriptor_type),
    values(NULL),
    nvalues(nvalues)
{
  // Copy values
  if (nvalues > 0) {
    this->values = new RNScalar [ nvalues ];
    for (int i = 0; i < nvalues; i++) {
      this->values[i] = values[i];
    }
  }
}



GSVDescriptor::
GSVDescriptor(const GSVDescriptor& descriptor)
  : descriptor_type(descriptor.descriptor_type),
    values(NULL),
    nvalues(descriptor.nvalues)
{
  // Copy values
  if (nvalues > 0) {
    values = new RNScalar [ nvalues ];
    for (int i = 0; i < nvalues; i++) {
      values[i] = descriptor.values[i];
    }
  }
}



GSVDescriptor::
~GSVDescriptor(void)
{
  // Delete values
  if (values) delete [] values;
}


RNScalar GSVDescriptor::
SquaredDistance(const GSVDescriptor& descriptor) const
{
  // Check descriptor type
  if (descriptor_type != descriptor.descriptor_type) return FLT_MAX;

  // Compute squared distance
  RNScalar sum = 0;
  for (int i = 0; i < nvalues; i++) {
    RNScalar delta = values[i] - descriptor.values[i];
    sum += delta * delta;
  }

  // Return sum of squared differences
  return sum;
}



static int
ComputeSpinImageDescriptor(GSVDescriptor& descriptor,
  const R2Grid& dh_position_x_grid, const R2Grid& dh_position_y_grid, int ix, int iy, 
  double height = 4, double radius = 2, int nstacks = 4, int nshells = 2)
{
  // Initialize descriptor
  if (descriptor.values) delete [] descriptor.values;
  descriptor.descriptor_type = GSV_NULL_DESCRIPTOR_TYPE;
  descriptor.values = NULL;
  descriptor.nvalues = 0;

  // Get spin image center position
  RNScalar cx = dh_position_x_grid.GridValue(ix, iy);
  if (cx == R2_GRID_UNKNOWN_VALUE) return 0;
  RNScalar cy = dh_position_y_grid.GridValue(ix, iy);
  if (cy == R2_GRID_UNKNOWN_VALUE) return 0;

  // Compute convenient variables
  int ixradius = (int) (radius * dh_position_x_grid.WorldToGridScaleFactor() + 0.5);
  int iyradius = (int) (height * dh_position_x_grid.WorldToGridScaleFactor() + 0.5);
  if (iyradius >= dh_position_y_grid.YResolution()) iyradius = dh_position_y_grid.YResolution();
  double radius_squared = radius * radius;

  // Allocate temporary values
  int nvalues = nshells * nstacks;
  RNScalar *values = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) values[i] = 0;

  // Count number of points in every spin image bin
  int count = 0;
  for (int i = -ixradius; i <= ixradius; i++) {
    if (ix + i < 0) continue;
    if (ix + i >= dh_position_x_grid.XResolution()) continue;
    for (int j = 0; j < iyradius; j++) {
      RNScalar x = dh_position_x_grid.GridValue(ix + i, j);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar y = dh_position_y_grid.GridValue(ix + i, j);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar dx = x - cx;
      RNScalar dy = y - cy;
      RNScalar rr = dx*dx + dy*dy;
      if (rr >= radius_squared) continue;
      RNScalar r = sqrt(rr);
      if (r >= radius) continue;
      int stack = nstacks * j / iyradius;
      int shell = (int) (nshells * (r / radius));
      values[stack*nshells + shell] += 1;
      count++;
    }
  }

  // Normalize number of points in spin image bins
  if (count > 0) {
    for (int i = 0; i < nvalues; i++) {
      values[i] /= count;
    }
  }

  // Update descriptor
  descriptor.descriptor_type = GSV_SPIN_IMAGE_DESCRIPTOR_TYPE;
  descriptor.nvalues = nvalues;
  descriptor.values = values;

  // Return success
  return 1;
}



static int
ComputeShapeContextDescriptor(GSVDescriptor& descriptor, const R2Image *image, 
  int cx, int cy, int radius = 16, int nsectors = 4, int nshells = 2)
{
  // Initialize descriptor
  if (descriptor.values) delete [] descriptor.values;
  descriptor.descriptor_type = GSV_NULL_DESCRIPTOR_TYPE;
  descriptor.values = NULL;
  descriptor.nvalues = 0;

  // Allocate temporary values
  int nvalues = nsectors * nshells;
  if (nvalues == 0) return 0;
  RNScalar *values = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) values[i] = 0;

  // Count number of points in every shape context bin
  double sum = 0;
  int radius_squared = radius * radius;
  for (int dx = -radius; dx <= radius; dx++) {
    int ix = cx + dx;
    if ((ix < 0) || (ix >= image->Width())) continue;
    for (int dy = -radius; dy <- radius; dy++) {
      int iy = cy + dy;
      if ((iy < 0) || (iy >= image->Width())) continue;
      RNScalar rr = dx*dx + dy*dy;
      if (rr >= radius_squared) continue;
      RNAngle angle = atan2(dy, dx) + RN_PI;
      int sector = (int) (nsectors * angle / RN_TWO_PI);
      RNScalar r = sqrt(rr);
      if (r >= radius) continue;
      int shell = (int) (nshells * (r / radius));
      if (sector >= nsectors) sector = nsectors-1;
      if (shell >= nshells) shell = nshells-1;
      RNRgb pixel = image->PixelRGB(ix, iy);
      RNScalar weight = pixel.Luminance();
      values[sector*nsectors + shell] += weight;
      sum += weight;
    }
  }

  // Normalize shape context bins
  if (sum > 0) {
    for (int i = 0; i < nvalues; i++) {
      values[i] /= sum;
    }
  }

  // Update descriptor
  descriptor.descriptor_type = GSV_SHAPE_CONTEXT_DESCRIPTOR_TYPE;
  descriptor.nvalues = nvalues;
  descriptor.values = values;

  // Return success
  return 1;
}



/////////////////////////////////////////////////////////////////////////
// GSVFeature member functions
////////////////////////////////////////////////////////////////////////

GSVFeature::
GSVFeature(void)
  : feature_type(GSV_NULL_FEATURE_TYPE),
    image(NULL),
    scanline(NULL),
    scan_point_index(-1),
    scan_position(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_direction(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_normal(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_scale(RN_UNKNOWN),
    image_position(RN_UNKNOWN, RN_UNKNOWN),
    image_direction(RN_UNKNOWN, RN_UNKNOWN),
    image_scale(RN_UNKNOWN),
    image_t(RN_UNKNOWN),
    score(RN_UNKNOWN),
    descriptor(),
    pixel_index(-1),
    index(-1)
{
}



GSVFeature::
GSVFeature(int feature_type, 
  GSVScanline *scanline, int scan_point_index, 
  const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, 
  RNScalar scan_scale, RNScalar score)
  : feature_type(feature_type),
    image(NULL),
    scanline(scanline),
    scan_point_index(scan_point_index),
    scan_position(scan_position),
    scan_direction(scan_direction),
    scan_normal(scan_normal),
    scan_scale(scan_scale),
    image_position(RN_UNKNOWN, RN_UNKNOWN),
    image_direction(RN_UNKNOWN, RN_UNKNOWN),
    image_scale(RN_UNKNOWN),
    image_t(RN_UNKNOWN),
    score(score),
    descriptor(),
    pixel_index(-1),
    index(-1)
{
}



GSVFeature::
GSVFeature(int feature_type, GSVImage *image, 
  const R2Point& image_position, const R2Vector& image_direction, 
  RNScalar image_scale, RNScalar score)
  : feature_type(feature_type),
    image(image),
    scanline(NULL),
    scan_point_index(-1),
    scan_position(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_direction(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_normal(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_scale(RN_UNKNOWN),
    image_position(image_position),
    image_direction(image_direction),
    image_scale(image_scale),
    image_t(RN_UNKNOWN),
    score(score),
    descriptor(),
    pixel_index(-1),
    index(-1)
{
  // Useful variables
  const R2Point unknown_point2(RN_UNKNOWN, RN_UNKNOWN);
  const R2Vector unknown_vector2(RN_UNKNOWN, RN_UNKNOWN);
  const double default_image_t = 10;

  // Compute scan stuff from image stuff
  if (image) {
    // Compute scanline
    if (!scanline) {
      RNScalar timestamp = image->Timestamp();
      GSVSegment *segment = image->Segment();
      if (!segment) return;
      if (segment->NScans() == 0) return;
      GSVScan *scan = segment->Scan(0);
      scanline = scan->FindScanlineBeforeTimestamp(timestamp);
    }

    // Compute scan position
    if (image_position != unknown_point2) {
      R3Ray ray = image->RayThroughUndistortedPosition(image_position);
      this->scan_position = ray.Point(default_image_t);
      this->image_t = 0;
    }

    // Compute scan direction
    if (image_direction != unknown_vector2) {
      R3Vector image_plane_direction = R3zero_vector;
      image_plane_direction += image_direction.X() * image->Pose().Right();
      image_plane_direction += image_direction.Y() * image->Pose().Up();
      this->scan_direction = image_plane_direction;
    }

    // Compute scan scale
    if (image_scale >= 0) {
      R2Point image_position2 = image_position + image_scale * image_direction;
      R3Ray ray2 = image->RayThroughUndistortedPosition(image_position2);
      R3Point scan_position2 = ray2.Point(default_image_t);
      this->scan_scale = R3Distance(this->scan_position, scan_position2);
    }
  }
}



GSVFeature::
GSVFeature(int feature_type, GSVImage *image, GSVScanline *scanline, int scan_point_index, 
  const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, RNScalar scan_scale,
  const R2Point& image_position, const R2Vector& image_direction, RNScalar image_scale, RNScalar image_t, 
  RNScalar score)
  : feature_type(feature_type),
    image(image),
    scanline(scanline),
    scan_point_index(scan_point_index),
    scan_position(scan_position),
    scan_direction(scan_direction),
    scan_normal(scan_normal),
    scan_scale(scan_scale),
    image_position(image_position),
    image_direction(image_direction),
    image_scale(image_scale),
    image_t(image_t),
    score(score),
    descriptor(),
    pixel_index(-1),
    index(-1)
{
  // Useful variables
  const R3Point unknown_point3(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);
  const R3Vector unknown_vector3(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);
  const R2Point unknown_point2(RN_UNKNOWN, RN_UNKNOWN);
  const R2Vector unknown_vector2(RN_UNKNOWN, RN_UNKNOWN);
  const double default_image_t = 10;

  // Compute scan stuff from image stuff
  if (image) {
    // Compute scanline
    if (!scanline) {
      RNScalar timestamp = image->Timestamp();
      GSVSegment *segment = image->Segment();
      if (!segment) return;
      if (segment->NScans() == 0) return;
      GSVScan *scan = segment->Scan(0);
      this->scanline = scan->FindScanlineBeforeTimestamp(timestamp);
    }

    // Compute scan position
    if ((scan_position == unknown_point3) && (image_position != unknown_point2)) {
      R3Ray ray = image->RayThroughUndistortedPosition(image_position);
      this->scan_position = ray.Point(default_image_t);
      this->image_t = 0;
    }

    // Compute scan direction
    if ((scan_direction == unknown_vector3) && (image_direction != unknown_vector2)) {
      R3Vector image_plane_direction = R3zero_vector;
      image_plane_direction += image_direction.X() * image->Pose().Right();
      image_plane_direction += image_direction.Y() * image->Pose().Up();
      this->scan_direction = image_plane_direction;
    }

    // Compute scan scale
    if ((scan_scale == RN_UNKNOWN) && (image_scale >= 0) && 
        (image_position != unknown_point2) && (image_direction != unknown_vector2))  {
      R2Point image_position2 = image_position + image_scale * image_direction;
      R3Ray ray2 = image->RayThroughUndistortedPosition(image_position2);
      R3Point scan_position2 = ray2.Point(default_image_t);
      this->scan_scale = R3Distance(this->scan_position, scan_position2);
    }
  }
}



void GSVFeature::
Draw(void) const
{
  // Parameters
  double radius = 0.1;

  // Check feature type
  if (feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
    R3Sphere(scan_position, radius).Draw();
  }
  else if (feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
    glLineWidth(5);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position - scan_scale * scan_direction);
    R3LoadPoint(scan_position + scan_scale * scan_direction);
    glEnd();
    glLineWidth(1);
  }
  else if (feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
    const R3Vector& normal = scan_normal;
    int dim = normal.MinDimension();
    R3Vector axis1 = R3xyz_triad[dim] % normal;
    axis1.Normalize();
    R3Vector axis2 = axis1 % normal;
    axis2.Normalize();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    glEnd();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    glEnd();
    glLineWidth(5);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position);
    R3LoadPoint(scan_position + normal);
    glEnd();
    glLineWidth(1);
  }
  else if (feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
    // Get pose
    R2Point distorted_position = image->DistortedPosition(image_position);
    int column_index = (int) (distorted_position.X() + 0.5);
    if (column_index < 0) column_index = 0;
    if (column_index >= image->Width()) column_index = image->Width()-1;
    R3Point viewpoint = image->Pose(column_index).Viewpoint();

    // Draw line
    glLineWidth(1);
    glBegin(GL_LINES);
    R3LoadPoint(viewpoint);
    R3LoadPoint(scan_position);
    glEnd();

    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);
  }
  else if (feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);
 
    // Draw line
    glLineWidth(3);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position);
    R3LoadPoint(scan_position + scan_scale*scan_direction);
    glEnd();
    glLineWidth(1);
  }
}



////////////////////////////////////////////////////////////////////////
// GSVFeatureSet member functions
////////////////////////////////////////////////////////////////////////

GSVFeatureSet::
GSVFeatureSet(GSVScene *scene)
  : scene(scene),
    features()
{
}


GSVFeatureSet::
GSVFeatureSet(const GSVFeatureSet& set)
  : scene(set.scene),
    features(set.features)
{
}


GSVFeatureSet& GSVFeatureSet::
operator=(const GSVFeatureSet& set)
{
  // Copy stuff
  scene = set.scene;
  features = set.features;
  return *this;
}



void GSVFeatureSet::
InsertFeature(GSVFeature *feature)
{
  // Insert feature
  features.Insert(feature);
}


void GSVFeatureSet::
RemoveFeature(GSVFeature *feature)
{
  // Remove feature
  features.Remove(feature);
}



void GSVFeatureSet::
Empty(void)
{
  // Remove all features
  features.Empty();
}



void GSVFeatureSet::
Draw(void) const
{
  // Draw features
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    feature->Draw();
  }
}



////////////////////////////////////////////////////////////////////////
// GSVCoprrespondence member functions
////////////////////////////////////////////////////////////////////////

GSVFeatureCorrespondence::
GSVFeatureCorrespondence(void)
  : score(RN_UNKNOWN),
    index(-1)
{
  // Assign features
  features[0] = NULL;
  features[1] = NULL;
}



GSVFeatureCorrespondence::
GSVFeatureCorrespondence(GSVFeature *feature0, GSVFeature *feature1, RNScalar score)
  : score(score),
    index(-1)
{
  // Assign features
  features[0] = feature0;
  features[1] = feature1;
}



void GSVFeatureCorrespondence::
Draw(void) const
{
  // Check features
  if (!features[0] || !features[1]) return;

  // Draw features
  features[0]->Draw();
  features[1]->Draw();
}



static int
GSVCompareCorrespondences(const void *data1, const void *data2)
{
  GSVFeatureCorrespondence *correspondence1 = *((GSVFeatureCorrespondence **) data1);
  GSVFeatureCorrespondence *correspondence2 = *((GSVFeatureCorrespondence **) data2);
  if (correspondence1->score > correspondence2->score) return -1;
  else if (correspondence1->score < correspondence2->score) return -1;
  else return 0;
}



////////////////////////////////////////////////////////////////////////
// GSVCoprrespondenceSet member functions
////////////////////////////////////////////////////////////////////////

GSVFeatureCorrespondenceSet::
GSVFeatureCorrespondenceSet(GSVScene *scene)
  : scene(scene),
    correspondences(),
    score(0)
{
}


GSVFeatureCorrespondenceSet::
GSVFeatureCorrespondenceSet(const GSVFeatureCorrespondenceSet& set)
  : scene(set.scene),
    correspondences(set.correspondences),
    score(set.score)
{
}


GSVFeatureCorrespondenceSet& GSVFeatureCorrespondenceSet::
operator=(const GSVFeatureCorrespondenceSet& set)
{
  // Copy everything
  scene = set.scene;
  correspondences = set.correspondences;
  score = set.score;
  return *this;
}



void GSVFeatureCorrespondenceSet::
InsertCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Insert correspondence
  correspondences.Insert(correspondence);
}



void GSVFeatureCorrespondenceSet::
RemoveCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Remove correspondence
  correspondences.Remove(correspondence);
}



void GSVFeatureCorrespondenceSet::
Empty(void)
{
  // Remove all correspondences
  correspondences.Empty();

  // Reset score
  score = 0;
}



void GSVFeatureCorrespondenceSet::
Sort(void)
{
  // Sort correspondences by score
  correspondences.Sort(GSVCompareCorrespondences);
}



void GSVFeatureCorrespondenceSet::
Draw(void) const
{
  // Draw correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    correspondence->Draw();
  }
}



////////////////////////////////////////////////////////////////////////
// GSVPath member functions
////////////////////////////////////////////////////////////////////////

GSVPath::
GSVPath(GSVSegment *segment, RNLength max_vertex_spacing)
  : segment(segment),
    spline(NULL),
    vertices(),
    index(-1)
{
  // Get scan
  if (segment->NScans() == 0) return;
  GSVScan *scan0 = segment->Scan(0);
  if (scan0->NScanlines() == 0) return;
  
  // Count total travel distance in scan
  RNScalar total_travel_distance = 0;
  R3Point previous_viewpoint = scan0->Scanline(0)->Pose().Viewpoint();
  for (int ie = 0; ie < scan0->NScanlines(); ie++) {
    GSVScanline *scanline = scan0->Scanline(ie);
    const R3Point& viewpoint = scanline->Pose().Viewpoint();
    total_travel_distance += R3Distance(viewpoint, previous_viewpoint);
    previous_viewpoint = viewpoint;
  }      
  
  // Compute spline control point parameters
  int num_vertices = 2 + (int) (total_travel_distance / max_vertex_spacing);
  if (num_vertices > scan0->NScanlines()) num_vertices = scan0->NScanlines();
  RNScalar vertex_spacing = total_travel_distance / (num_vertices - 1);
  
  // Create spline vertices
  RNScalar travel_distance = 0;
  previous_viewpoint = scan0->Scanline(0)->Pose().Viewpoint();
  int previous_vertex_scanline_index = 0;
  RNScalar previous_vertex_travel_distance = 0;
  RNScalar next_vertex_travel_distance = vertex_spacing;
  RNArray<R3Point *> vertex_positions;
  for (int ie = 0; ie < scan0->NScanlines(); ie++) {
    GSVScanline *scanline = scan0->Scanline(ie);
    const R3Point& viewpoint = scanline->Pose().Viewpoint();
    travel_distance += R3Distance(viewpoint, previous_viewpoint);
    previous_viewpoint = viewpoint;
    if ((ie == 0) || (ie == scan0->NScanlines()-1) || 
        ((travel_distance >= next_vertex_travel_distance) && 
         (total_travel_distance - travel_distance > 1.5 * vertex_spacing))) {
      // Add spline vertex
      GSVPathVertex *vertex = new GSVPathVertex();
      vertex->path = this;
      vertex->pose = scanline->Pose();
      vertex->translation = R3zero_vector;
      vertex->rotation = R3zero_vector;
      vertex->parameter = vertex_positions.NEntries();
      vertex->timestamp = scanline->Timestamp();
      vertices.Insert(vertex);
      vertex_positions.Insert((R3Point *) &(scanline->Pose().Viewpoint()));
      next_vertex_travel_distance = travel_distance + vertex_spacing;
      previous_vertex_travel_distance = travel_distance;
      previous_vertex_scanline_index = ie;
    }
  }

  // Create spline
  spline = new R3CatmullRomSpline(vertex_positions);
}



GSVPath::
~GSVPath(void)
{
  // Delete spline
  if (spline) delete spline;
}



GSVPose GSVPath::
TransformedPose(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return GSVnull_pose;
  if (NVertices() == 1) return VertexPose(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector translation0 = VertexTranslation(iu0);
  R3Vector translation1 = VertexTranslation(iu1);
  R3Vector rotation0 = VertexRotation(iu0);
  R3Vector rotation1 = VertexRotation(iu1);
  R3Vector translation = (1-t) * translation0 + t * translation1;
  R3Vector rotation = (1-t) * rotation0 + t * rotation1;

  // Return transformed pose at parameter value
  GSVPose pose = Pose(u);
  pose.Translate(translation);
  pose.Rotate(rotation);
  return pose;
}



GSVPose GSVPath::
Pose(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return GSVnull_pose;
  if (NVertices() == 1) return VertexPose(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  GSVPose pose0 = VertexPose(iu0);
  GSVPose pose1 = VertexPose(iu1);
  return GSVInterpolatedPose(pose0, pose1, t);
}



R3Vector GSVPath::
Translation(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return R3zero_vector;
  if (NVertices() == 1) return VertexTranslation(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector translation0 = VertexTranslation(iu0);
  R3Vector translation1 = VertexTranslation(iu1);
  return (1-t) * translation0 + t * translation1;
}



R3Vector GSVPath::
Rotation(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return R3zero_vector;
  if (NVertices() == 1) return VertexRotation(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector rotation0 = VertexRotation(iu0);
  R3Vector rotation1 = VertexRotation(iu1);
  return (1-t) * rotation0 + t * rotation1;
}



RNScalar GSVPath::
Timestamp(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return 0.0;
  if (NVertices() == 1) return VertexTimestamp(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  RNScalar timestamp0 = VertexTimestamp(iu0);
  RNScalar timestamp1 = VertexTimestamp(iu1);
  return (1-t) * timestamp0 + t * timestamp1;
}



RNScalar GSVPath::
Parameter(RNScalar timestamp) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return RNScalar(0.0);
  if (NVertices() == 1) return VertexParameter(0);

  // Return parameter u at timestamp
  int i0 = FindVertexIndexBeforeTimestamp(timestamp);
  assert((i0 >= 0) && (i0 < NVertices()-1));
  RNScalar timestamp0 = VertexTimestamp(i0);
  RNScalar timestamp1 = VertexTimestamp(i0 + 1);
  RNScalar denom = timestamp1 - timestamp0;
  if (RNIsZero(denom)) return 0.5 * (timestamp0 + timestamp1);
  RNScalar t = (timestamp - timestamp0) / denom;
  RNScalar u0 = VertexParameter(i0);
  RNScalar u1 = VertexParameter(i0 + 1);
  return RNScalar((1-t)*u0 + t*u1);
}



RNInterval GSVPath::
ParameterRange(void) const
{
  // Return range of parameter u
  GSVPathVertex *tail = vertices.Tail();
  return RNInterval(0.0, tail->parameter);
}



RNInterval GSVPath::
TimestampRange(void) const
{
  // Return range of timestamps
  GSVPathVertex *tail = vertices.Tail();
  return RNInterval(0.0, tail->timestamp);
}



int GSVPath::
NVertices(void) const
{
  // Return number of vertices
  return vertices.NEntries();
}



GSVPathVertex *GSVPath::
Vertex(int k) const
{
  // Return kth vertex
  return vertices.Kth(k);
}



const GSVPose& GSVPath::
VertexPose(int vertex_index) const
{
  // Return pose at kth vertex
  return vertices[vertex_index]->pose;
}



const R3Point& GSVPath::
VertexPosition(int vertex_index) const
{
  // Return viewpoint position at kth vertex
  return vertices[vertex_index]->pose.Viewpoint();
}



const R3Vector& GSVPath::
VertexTranslation(int vertex_index) const
{
  // Return translation at kth vertex
  return vertices[vertex_index]->translation;
}



const R3Vector& GSVPath::
VertexRotation(int vertex_index) const
{
  // Return rotation at kth vertex
  return vertices[vertex_index]->rotation;
}



RNScalar GSVPath::
VertexParameter(int vertex_index) const
{
  // Return parameter u at kth vertex
  return vertices[vertex_index]->parameter;
}



RNScalar GSVPath::
VertexTimestamp(int vertex_index) const
{
  // Return timestamp at kth vertex
  return vertices[vertex_index]->timestamp;
}



GSVPose GSVPath::
VertexTransformedPose(int vertex_index) const
{
  // Return transformed pose 
  R3Vector translation = VertexTranslation(vertex_index);
  R3Vector rotation = VertexRotation(vertex_index);
  GSVPose pose = VertexPose(vertex_index);
  pose.Translate(translation);
  pose.Rotate(rotation);
  return pose;
}



void GSVPath::
Draw(void) const
{
  // Draw spline
  if (!spline) return;
  spline->Draw();
}


int GSVPath::
FindVertexIndexBeforeTimestamp(RNScalar timestamp) const
{
  // Binary search
  if (NVertices() == 0) return -1;
  if (timestamp <= VertexTimestamp(0)) return 0;
  if (timestamp >= VertexTimestamp(NVertices()-2)) return NVertices()-2;
  return FindVertexIndexBeforeTimestamp(timestamp, 0, NVertices()-2);
}



int GSVPath::
FindVertexIndexBeforeTimestamp(RNScalar timestamp, int imin, int imax) const
{
  // Binary search
  int i = (imin + imax) / 2;
  if (i == imin) return imin;
  assert((i >= 0) && (i < imax));
  RNScalar t = VertexTimestamp(i);
  if (t > timestamp) return FindVertexIndexBeforeTimestamp(timestamp, imin, i);
  else if (t < timestamp) return FindVertexIndexBeforeTimestamp(timestamp, i, imax);
  else return i;
}



GSVPathVertex:: 
GSVPathVertex(void)
  : path(NULL),
    pose(0,0,0,0,0,0,0),
    translation(0,0,0),
    rotation(0,0,0),
    parameter(-1),
    timestamp(-1),
    index(-1)
{
}



////////////////////////////////////////////////////////////////////////
// GSVPoseOptimization member functions
////////////////////////////////////////////////////////////////////////

GSVPoseOptimization::
GSVPoseOptimization(GSVScene *scene)
  : scene(scene),
    features(scene),
    correspondences(scene),
    lasers(),
    cameras(),
    segments(),
    scanlines(),
    images(),
    vertices(),
    score(-1),
    pixels(),
    laser_nv(0),
    camera_nv(0),
    vertex_nv(0),
    pixel_nv(0)
{
  // Set feature extraction parameters
  sift_image_scale = 0.125;

  // Set correspondence parameters
  create_point_point_correspondences = TRUE;
  create_plane_plane_correspondences = TRUE;
  create_pixel_point_correspondences = TRUE;
  create_pixel_plane_correspondences = TRUE;
  create_pixel_pixel_correspondences = TRUE;
  min_intracorrespondence_path_parameter_difference = 1;
  max_intracorrespondence_euclidean_distance = 1;
  max_intracorrespondence_spin_image_descriptor_distance = 0.25;
  max_intracorrespondence_sift_descriptor_distance = 0;
  max_intracorrespondence_line_descriptor_distance = 0;
  max_intracorrespondence_direction_angle = RN_PI / 6.0;
  max_intracorrespondence_normal_angle = RN_PI / 6.0;
  min_intrasegment_path_parameter_difference = 10;
  max_intrasegment_distance_ratio = 0.1;
  icp_max_distance_start = 20;
  icp_max_distance_end = 1.0;

  // Set optimization parameters
  RNScalar laser_translation_inertia = 1E6;
  RNScalar camera_translation_inertia = 1E6;
  RNScalar vertex_translation_inertia = 1E-3;
  laser_inertia_weights[TX] = laser_translation_inertia;
  laser_inertia_weights[TY] = laser_translation_inertia;
  laser_inertia_weights[TZ] = laser_translation_inertia;
  camera_inertia_weights[TX] = camera_translation_inertia;
  camera_inertia_weights[TY] = camera_translation_inertia;
  camera_inertia_weights[TZ] = camera_translation_inertia;
  vertex_inertia_weights[TX] = vertex_translation_inertia;
  vertex_inertia_weights[TY] = vertex_translation_inertia;
  vertex_inertia_weights[TZ] = vertex_translation_inertia;
  vertex_inertia_weights[RX] = 1E6 * vertex_translation_inertia;  // don't tilt
  vertex_inertia_weights[RY] = 1E6 * vertex_translation_inertia;  // don't tilt
  vertex_inertia_weights[RZ] = 1E2 * vertex_translation_inertia; // 100 meters ~ 1 radian ~ 60 degrees
  pixel_inertia_weights[T] = 1E-9;
  intra_segment_rigidity_weight = 1;
  intra_segment_rigidity_radius = 4;
  intra_segment_rigidity_sigma = 0.5 * intra_segment_rigidity_radius;
  intra_segment_rigidity_sigma_squared = intra_segment_rigidity_sigma * intra_segment_rigidity_sigma;
  scan_point_scan_point_correspondence_weight = (2*intra_segment_rigidity_radius + 1) * intra_segment_rigidity_weight;
  scan_point_scan_plane_correspondence_weight = scan_point_scan_point_correspondence_weight;
  scan_point_image_point_correspondence_weight = scan_point_scan_point_correspondence_weight;
  image_point_scan_plane_correspondence_weight = scan_point_image_point_correspondence_weight;
  image_point_image_point_correspondence_weight = scan_point_image_point_correspondence_weight;
  max_vertex_spacing = 2;

  // Initialize optimization variable selection
  for (int i = 0; i < LASER_NV; i++) laser_v[i] = -1;
  for (int i = 0; i < CAMERA_NV; i++) camera_v[i] = -1;
  for (int i = 0; i < VERTEX_NV; i++) vertex_v[i] = -1;
  for (int i = 0; i < PIXEL_NV; i++) pixel_v[i] = -1;

  // Check scene
  if (scene) {
    // Create data
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);

      // Create data for every laser
      for (int is = 0; is < run->NLasers(); is++) {
        GSVLaser *laser = run->Laser(is);
        assert(!laser->Data());
        LaserData *laser_data = new LaserData();
        laser_data->laser = laser;
        laser_data->translation = R3zero_vector;
        laser_data->index = lasers.NEntries();
        laser->SetData(laser_data);
        lasers.Insert(laser);
      }

      // Create data for every camera
      for (int is = 0; is < run->NCameras(); is++) {
        GSVCamera *camera = run->Camera(is);
        assert(!camera->Data());
        CameraData *camera_data = new CameraData();
        camera_data->camera = camera;
        camera_data->translation = R3zero_vector;
        camera_data->index = cameras.NEntries();
        camera->SetData(camera_data);
        cameras.Insert(camera);
      }

      // Create data for every segment
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        
        // Create path
        GSVPath *path = new GSVPath(segment, max_vertex_spacing);
        if (!path || (path->NVertices() == 0)) {
          fprintf(stderr, "Unable to create path for segment %d\n", is);
          continue;
        }

        // Insert vertices from path
        for (int i = 0; i < path->NVertices(); i++) {
          GSVPathVertex *vertex = path->Vertex(i);
          vertex->index = vertices.NEntries();
          vertices.Insert(vertex);
        }

        // Create segment data
        assert(!segment->Data());
        SegmentData *segment_data = new SegmentData();
        segment_data->segment = segment;
        segment_data->path = path;
        segment_data->index = segments.NEntries();
        segment->SetData((void *) segment_data);
        segments.Insert(segment);

        // Create scanline data 
        for (int ia = 0; ia < segment->NScans(); ia++) {
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            R3Point viewpoint = scanline->Pose().Viewpoint();
            RNScalar timestamp = scanline->Timestamp();
            RNScalar path_parameter = path->Parameter(timestamp);
            R3Point path_viewpoint = path->Pose(path_parameter).Viewpoint();
            assert(!scanline->Data());
            ScanlineData *scanline_data = new ScanlineData();
            scanline_data->scanline = scanline;
            scanline_data->path_parameter = path_parameter;
            scanline_data->index = scanlines.NEntries();
            scanline->SetData((void *) scanline_data);
            scanlines.Insert(scanline);
          }
        }

        // Create image data
        for (int ip = 0; ip < segment->NPanoramas(); ip++) {
          GSVPanorama *panorama = segment->Panorama(ip);
          for (int ii = 0; ii < panorama->NImages(); ii++) {
            GSVImage *image = panorama->Image(ii);
            R3Point viewpoint = image->Pose().Viewpoint();
            RNScalar timestamp = image->Timestamp();
            RNScalar path_parameter = path->Parameter(timestamp);
            R3Point path_viewpoint = path->Pose(path_parameter).Viewpoint();
            assert(!image->Data());
            ImageData *image_data = new ImageData();
            image_data->image = image;
            image_data->path_parameter = path_parameter;
            image_data->index = images.NEntries();
            image->SetData((void *) image_data);
            images.Insert(image);
          }
        }
      }
    }
  }
}



R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVScanline *scanline) const
{
  // Get convenient variables
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return R3identity_affine;
  GSVScan *scan = scanline->Scan();
  if (!scan) return R3identity_affine;
  GSVLaser *laser = scanline->Scan()->Laser();
  if (!laser) return R3identity_affine;
  LaserData *laser_data = (LaserData *) laser->Data();
  if (!laser_data) return R3identity_affine;
  GSVSegment *segment = scan->Segment();
  if (!segment) return R3identity_affine;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return R3identity_affine;
  GSVPath *path = segment_data->path;
  if (!path) return R3identity_affine;

  // Get path pose
  const GSVPose& path_pose = path->Pose(scanline_data->path_parameter);
  const R3Point& path_origin = path_pose.Viewpoint();
  const R3Vector path_right = path_pose.Right();
  const R3Vector path_towards = path_pose.Towards();
  const R3Vector path_up = path_pose.Up();
  const R3Vector& path_translation = path->Translation(scanline_data->path_parameter);
  const R3Vector& path_rotation = path->Rotation(scanline_data->path_parameter);

  // Convert laser transformation from path coordinate system to world coordinate system
  R3Vector laser_translation;
  R3Vector lt = laser_data->translation;
  laser_translation[0] = lt.X()*path_right.X() + lt.Y()*path_towards.X() + lt.Z()*path_up.X();
  laser_translation[1] = lt.X()*path_right.Y() + lt.Y()*path_towards.Y() + lt.Z()*path_up.Y();
  laser_translation[2] = lt.X()*path_right.Z() + lt.Y()*path_towards.Z() + lt.Z()*path_up.Z();

  // Return transformation 
  R4Matrix matrix(R4identity_matrix);
  matrix.Translate(path_translation);
  matrix.Translate(path_origin.Vector());
  matrix.Rotate(path_rotation);
  matrix.Translate(-(path_origin.Vector()));
  matrix.Translate(laser_translation);
  return R3Affine(matrix);
}



R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVImage *image) const
{
  // Get convenient variables
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return R3identity_affine;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return R3identity_affine;
  GSVCamera *camera = image->Tapestry()->Camera();
  if (!camera) return R3identity_affine;
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return R3identity_affine;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return R3identity_affine;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return R3identity_affine;
  GSVPath *path = segment_data->path;
  if (!path) return R3identity_affine;

  // Get path pose
  const GSVPose& path_pose = path->Pose(image_data->path_parameter);
  const R3Point& path_origin = path_pose.Viewpoint();
  const R3Vector path_right = path_pose.Right();
  const R3Vector path_towards = path_pose.Towards();
  const R3Vector path_up = path_pose.Up();
  const R3Vector& path_translation = path->Translation(image_data->path_parameter);
  const R3Vector& path_rotation = path->Rotation(image_data->path_parameter);

  // Convert laser transformation from path coordinate system to world coordinate system
  R3Vector camera_translation;
  R3Vector lt = camera_data->translation;
  camera_translation[0] = lt.X()*path_right.X() + lt.Y()*path_towards.X() + lt.Z()*path_up.X();
  camera_translation[1] = lt.X()*path_right.Y() + lt.Y()*path_towards.Y() + lt.Z()*path_up.Y();
  camera_translation[2] = lt.X()*path_right.Z() + lt.Y()*path_towards.Z() + lt.Z()*path_up.Z();

  // Return transformation 
  R4Matrix matrix(R4identity_matrix);
  matrix.Translate(path_translation);
  matrix.Translate(path_origin.Vector());
  matrix.Rotate(path_rotation);
  matrix.Translate(-(path_origin.Vector()));
  matrix.Translate(camera_translation);
  return R3Affine(matrix);
}



R3Point GSVPoseOptimization::
FeaturePosition(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature position (possibly with optimization transformation applied)
  R3Point result = feature->scan_position;
  if (!apply_pose_transformation) return result;
  else if (feature->image) result.Transform(OptimizedTransformation(feature->image));
  else if (feature->scanline) result.Transform(OptimizedTransformation(feature->scanline));
  return result;
}



R3Vector GSVPoseOptimization::
FeatureDirection(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature direction (possibly with optimization transformation applied)
  R3Vector result = feature->scan_direction;
  if (!apply_pose_transformation) return result;
  else if (feature->image) result.Transform(OptimizedTransformation(feature->image));
  else if (feature->scanline) result.Transform(OptimizedTransformation(feature->scanline));
  return result;
}



R3Vector GSVPoseOptimization::
FeatureNormal(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature normal (possibly with optimization transformation applied)
  R3Vector result = feature->scan_normal;
  if (!apply_pose_transformation) return result;
  else if (feature->image) result.Transform(OptimizedTransformation(feature->image));
  else if (feature->scanline) result.Transform(OptimizedTransformation(feature->scanline));
  return result;
}



RNRgb GSVPoseOptimization::
FeatureColor(GSVFeature *feature) const
{
  // Colors
  const int max_feature_colors = 6;
  const RNRgb feature_colors[max_feature_colors] = {
    RNRgb(1, 0, 0), RNRgb(0, 1, 0), RNRgb(0, 0, 1), 
    RNRgb(0.8, 0.8, 0), RNRgb(0, 0.8, 0.8), RNRgb(0.8, 0, 0.8)
  };

  // Find run index and segment index
  GSVScanline *scanline = feature->scanline;
  if (!scanline) return RNblack_rgb;
  GSVScan *scan = scanline->Scan();
  if (!scan) return RNblack_rgb;
  GSVSegment *segment = scan->Segment();
  if (!segment) return RNblack_rgb;
  GSVRun *run = segment->Run();
  if (!run) return RNblack_rgb;
  int segment_index = segment->RunIndex();
  int run_index = run->SceneIndex();

  // Return default feature color
  return feature_colors[(7*run_index + segment_index) % max_feature_colors];
}



void GSVPoseOptimization::
InsertFeature(GSVFeature *feature)
{
  // Insert into list of features
  assert(feature->index == -1);
  feature->index = features.features.NEntries();
  features.features.Insert(feature);
}


void GSVPoseOptimization::
RemoveFeature(GSVFeature *feature)
{
  // Remove from list of features
  assert(feature->index >= 0);
  RNArrayEntry *entry = features.features.KthEntry(feature->index);
  GSVFeature *tail = features.features.Tail();
  tail->index = feature->index;
  features.features.EntryContents(entry) = tail;
  features.features.RemoveTail();
  feature->index = -1;
}



void GSVPoseOptimization::
InsertCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Insert correspondence
  assert(correspondence->index == -1);
  correspondence->index = correspondences.correspondences.NEntries();
  correspondences.correspondences.Insert(correspondence);
}


void GSVPoseOptimization::
RemoveCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Remove correspondence
  RNArrayEntry *entry = correspondences.correspondences.KthEntry(correspondence->index);
  GSVFeatureCorrespondence *tail = correspondences.correspondences.Tail();
  tail->index = correspondence->index;
  correspondences.correspondences.EntryContents(entry) = tail;
  correspondences.correspondences.RemoveTail();
  correspondence->index = -1;
}



////////////////////////////////////////////////////////////////////////
// Harris corner feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CreateHarrisCornerFeatures(GSVPoseOptimization *optimization, GSVImage *image)
{
  // Parameters
  int radius = 8;
  RNScalar kappa = 0.05;
  RNScalar threshold = 0.001;
  int min_spacing = 16;

  // Get image
  R2Image *img = image->UndistortedImage();
  if (!img) return 0;

  // Create grid with corners 
  R2Grid grid(*img, 0);
  grid.HarrisCornerFilter(radius, kappa);
  grid.Threshold(threshold, 0, R2_GRID_KEEP_VALUE);
  grid.MaskNonMaxima(min_spacing);

  // Create features
  for (int i = 0; i < grid.XResolution(); i++) {
    for (int j = 0; j < grid.YResolution(); j++) {
      // Get info
      R2Vector direction(0, 0);
      R2Point position(i, j);
      RNScalar score = grid.GridValue(i, j);
      if ((score == 0) || (score == R2_GRID_UNKNOWN_VALUE)) continue;
      GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, image, position, direction, radius, score);
      if (!feature) { fprintf(stderr, "Unable to create harris corner feature\n"); return 0; }
      ComputeShapeContextDescriptor(feature->descriptor, img, i, j);
      optimization->InsertFeature(feature);
    }
  }

#if 0
  char buffer[128];
  R2Grid tmp(grid);
  tmp.Dilate(2);
  tmp.Multiply(10);
  tmp.Add(R2Grid(*img, 0));
  sprintf(buffer, "mkdir -p tmp");
  system(buffer);
  sprintf(buffer, "tmp/%d.grd", image->SceneIndex());
  tmp.Write(buffer);
  printf("%d %d\n", image->SceneIndex(), grid->Cardinality());
#endif

  // Delete image
  if (img) delete img;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageCornerFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/image_corner_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Make run directories
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
  }

  // Create corner features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          // int saved_nfeatures = NFeatures();
          if (!CreateHarrisCornerFeatures(this, image)) return 0;
          // printf("%d %d %d %d : %d\n", ir, is, ip, ii, NFeatures() - saved_nfeatures);
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Sift feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CreateSiftFeatures(GSVPoseOptimization *optimization, GSVImage *image, const char *directory_name)
{
  // Get filenames
  char pgm_name[4096], sift_name[4096];
  int image_index = image->PanoramaIndex();
  GSVPanorama *panorama = image->Panorama();
  int panorama_index = panorama->SegmentIndex();
  GSVSegment *segment = panorama->Segment();
  int segment_index = segment->RunIndex();
  sprintf(pgm_name, "%s/%02d_%06d_%02d.pgm", directory_name, segment_index, panorama_index, image_index);
  sprintf(sift_name, "%s/%02d_%06d_%02d.sft", directory_name, segment_index, panorama_index, image_index);

  // Check if sift file already exists 
  if (!RNFileExists(sift_name)) {
    // Get image
    R2Image *img = image->UndistortedImage();
    if (!img) return 0;
    
    // Scale image (this could be better)
    RNScalar scale = optimization->sift_image_scale;
    int width = scale * img->Width();
    int height = scale * img->Height();
    R2Image scaled_img(width, height, 3);
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        int x = (int) (i / scale + 0.5);
        int y = (int) (j / scale + 0.5);
        if ((x < 0) || (x >= img->Width())) continue;
        if ((y < 0) || (y >= img->Height())) continue;
        RNRgb pixel = img->PixelRGB(x, y);
        scaled_img.SetPixelRGB(i, j, pixel);
      }
    }
    

    // Write pgm image
    if (!scaled_img.Write(pgm_name)) {
      delete img;
      return 0;
    }

    // Run program to extract sift features from image
    char sift_command[4096];
    sprintf(sift_command, "sift < %s > %s", pgm_name, sift_name);
    system(sift_command);
    sprintf(sift_command, "rm -f %s", pgm_name);
    system(sift_command);
    
    // Delete image
    delete img;
  }

  // Read sift file
  //return optimization->ReadSiftFile(image, sift_name);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageSiftFeatures(void)
{
  // Make temporary directory
  //  char tmp_directory[4096], mkdir_cmd[4096];
  //sprintf(tmp_directory, "%s/image_sift_features", scene->CacheDataDirectoryName());
  //sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  //system(mkdir_cmd);

  // Make run directories
  /*for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
    }*/

  // Create sift features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    //char run_directory[4096];
    //sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
	  R2Image *img = image->UndistortedImage();
	  delete img;
	  //if (!CreateSiftFeatures(this, image, run_directory)) return 0;
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Line feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CreateLineFeatures(GSVPoseOptimization *optimization, GSVImage *image, const char *directory_name)
{
  // Get filenames
  char image_name[4096], line_name[4096];
  int image_index = image->PanoramaIndex();
  GSVPanorama *panorama = image->Panorama();
  int panorama_index = panorama->SegmentIndex();
  GSVSegment *segment = panorama->Segment();
  int segment_index = segment->RunIndex();
  sprintf(image_name, "%s/%02d_%06d_%02d.pgm", directory_name, segment_index, panorama_index, image_index);
  sprintf(line_name, "%s/%02d_%06d_%02d.lin", directory_name, segment_index, panorama_index, image_index);

  // TEMPORARY
  if (!RNFileExists(line_name)) return 1;

  // Check if line file already exists 
  if (!RNFileExists(line_name)) {
    // Get image
    R2Image *img = image->UndistortedImage();
    if (!img) return 0;

    // Write image to temporary directory
    if (!img->Write(image_name)) {
      delete img;
      return 0;
    }

    // Run program to extract line features from image
    char line_command[4096];
    sprintf(line_command, "pfm2lin %s %s", image_name, line_name);
    system(line_command);
    sprintf(line_command, "rm -f %s", image_name);
    system(line_command);
    
    // Delete image
    delete img;
  }

  // Read line file
  return optimization->ReadLineFile(image, line_name);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageLineFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/image_line_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Make run directories
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
  }

  // Create line features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          if (!CreateLineFeatures(this, image, run_directory)) return 0;
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Pole feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CompareScalarPointers(const void *data1, const void *data2)
{
  const RNScalar *ptr1 = *((const RNScalar **) data1);
  const RNScalar *ptr2 = *((const RNScalar **) data2);
  RNScalar value1 = *ptr1;
  RNScalar value2 = *ptr2;
  if (value1 > value2) return -1;
  else if (value2 > value1) return 1;
  else return 0;
}



static void 
FilterVotes(R2Grid& vote_grid, const R2Grid& depth_grid, 
  int min_iy, int max_iy, RNScalar y_sigma, RNScalar depth_sigma)
{
  // Make copy of vote grid
  R2Grid vote_copy(vote_grid);
  vote_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  
  // Get convenient variables
  const RNScalar sqrt_two_pi = sqrt(RN_TWO_PI);
  double y_fac = 1.0 / (sqrt_two_pi * y_sigma);
  double y_denom = -2.0 * y_sigma * y_sigma;
  double d_fac = 1.0 / (sqrt_two_pi * depth_sigma);
  double d_denom = -2.0 * depth_sigma * depth_sigma;
  int y_radius = (int) (3.0*y_sigma) + 1;
  RNScalar depth_radius = 3 * depth_sigma;

  // Set every vote to be bilateral filter of vertical region
  for (int ix = 0; ix < vote_grid.XResolution(); ix++) {
    for (int iy = min_iy; iy <= max_iy; iy++) {
      RNScalar vote = vote_copy.GridValue(ix, iy);
      if (vote == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar depth = depth_grid.GridValue(ix, iy);
      if (depth == R2_GRID_UNKNOWN_VALUE) continue;

      // Compute blurred vote value
      RNScalar sum = 0;
      RNScalar weight = 0;
      for (int y = iy - y_radius; y <= iy + y_radius; y++) {
        if ((y < 0) || (y >= vote_grid.YResolution())) continue;
        RNScalar vote_sample = vote_copy.GridValue(ix, y);
        if (vote_sample == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar depth_sample = depth_grid.GridValue(ix, y);
        if (depth_sample == R2_GRID_UNKNOWN_VALUE) continue;
        int dy = y - iy;
        RNScalar dd = depth - depth_sample;
        if (dd > depth_radius) continue;
        RNScalar w = y_fac * exp(dy*dy/y_denom) * d_fac * exp(dd*dd/d_denom);
        sum += w * vote_sample;
        weight += w;
      }
      
      // Set grid value
      // if (weight > 0) vote_grid.SetGridValue(ix, iy, sum / weight);
      vote_grid.SetGridValue(ix, iy, sum);
    }
  }
}



int GSVPoseOptimization::
CreateScanPoleFeatures(void)
{
  // Parameters
  const RNScalar min_height = 0.25;
  const RNScalar max_height = 4.0;
  const RNScalar max_depth = 10;
  const RNScalar min_depth_difference = 0.5;
  const RNScalar expected_pole_radius = 0.25;
  const RNScalar expected_pole_height = 2.0;
  const RNScalar min_feature_spacing = 1;
  const RNScalar min_coverage = 0.5;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create pole features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;

        // Read points (needed for EstimatedGroundZ)
        if (!scan->ReadPoints()) {
          fprintf(stderr, "Unable to read points for scan\n");
          return 0;
        }

        // Read DH grids
        char image_name[4096];
        R2Grid dh_viewpoint_depth_grid;
        R2Grid dh_scanline_index_grid;
        R2Grid dh_point_index_grid;
        R2Grid dh_position_x_grid;
        R2Grid dh_position_y_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_ViewpointDepth.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_viewpoint_depth_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_y_grid.Read(image_name)) return 0;

        // Create DH grid containing vertical valley indicator function
        int pole_radius = (int) (dh_viewpoint_depth_grid.WorldToGridScaleFactor() * expected_pole_radius) + 1;
        R2Grid dh_vertical_valley_grid(dh_viewpoint_depth_grid);
        dh_vertical_valley_grid.Clear(-1);
        for (int i = 0; i < dh_viewpoint_depth_grid.XResolution(); i++) {
          for (int j = 0; j < dh_viewpoint_depth_grid.YResolution(); j++) {
            int left_ix = i - pole_radius;
            int right_ix = i + pole_radius;
            if (left_ix < 0) left_ix = 0;
            if (right_ix > dh_viewpoint_depth_grid.XResolution()-1) 
              right_ix = dh_viewpoint_depth_grid.XResolution()-1;
            RNScalar value = dh_viewpoint_depth_grid.GridValue(i, j);
            if (value == R2_GRID_UNKNOWN_VALUE) continue;
            if (value > max_depth) continue;
            RNScalar left_value = dh_viewpoint_depth_grid.GridValue(left_ix, j);
            if (left_value == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar right_value = dh_viewpoint_depth_grid.GridValue(right_ix, j);
            if (right_value == R2_GRID_UNKNOWN_VALUE) continue;
            if (left_value - value < min_depth_difference) continue;
            if (right_value - value < min_depth_difference) continue;
            dh_vertical_valley_grid.SetGridValue(i, j, 1);
          }
        }

        // Create DH vote grid
        R2Grid dh_vote_grid(dh_vertical_valley_grid);
        int y1 = (int) (dh_vote_grid.WorldToGridScaleFactor() * min_height);
        int y2 = (int) (dh_vote_grid.WorldToGridScaleFactor() * max_height);
        if (y1 >= dh_vote_grid.YResolution()) y1 = dh_vote_grid.YResolution()-1;
        if (y2 >= dh_vote_grid.YResolution()) y2 = dh_vote_grid.YResolution()-1;
        double y_sigma = 0.5 * dh_vote_grid.WorldToGridScaleFactor() * expected_pole_height;
        double depth_sigma = expected_pole_radius;
        FilterVotes(dh_vote_grid, dh_viewpoint_depth_grid, y1, y2, y_sigma, depth_sigma);

        // Find local maxima of dh_vote_grid
        RNArray<const RNScalar *> local_maxima;
        R2Grid dh_maxima_grid(dh_vote_grid);
        dh_maxima_grid.MaskNonMaxima(2 * expected_pole_radius * dh_maxima_grid.WorldToGridScaleFactor());
        dh_maxima_grid.Threshold(min_coverage, R2_GRID_UNKNOWN_VALUE, R2_GRID_KEEP_VALUE);
        for (int i = 0; i < dh_maxima_grid.XResolution(); i++) {
          // Find maximum value in column
          int maximum_j = -1;
          RNScalar maximum_value = -FLT_MAX;
          for (int j = y1; j <= y2; j++) {
            RNScalar value = dh_maxima_grid.GridValue(i, j);
            if (value == R2_GRID_UNKNOWN_VALUE) continue;
            if (value == 0) continue;
            if (value < maximum_value) continue;
            maximum_value = value;
            maximum_j = j;
          }

          // Insert local maximum
          if (maximum_j >= 0) {
            const RNScalar *local_maximum = &dh_maxima_grid(i, maximum_j);
            local_maxima.Insert(local_maximum);
          }
        }

        // Sort local maxima
        local_maxima.Sort(CompareScalarPointers);

        // Create features at local maxima with minimum spacing
        RNArray<GSVFeature *> created_features;
        RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;
        for (int k = 0; k < local_maxima.NEntries(); k++) {
          // Get local maximum info
          const RNScalar *local_maximum = local_maxima.Kth(k);
          int grid_index = local_maximum - dh_maxima_grid.GridValues();
          assert((grid_index >= 0) && (grid_index < dh_maxima_grid.NEntries()));

          // Determine the x,y position
          RNScalar x = dh_position_x_grid.GridValue(grid_index);
          RNScalar y = dh_position_y_grid.GridValue(grid_index);

          // Check if there is a previously created (stronger) feature within minimum spacing
          RNBoolean well_spaced = TRUE;
          R2Point local_maximum_position(x, y);
          for (int j = 0; j < created_features.NEntries(); j++) {
            GSVFeature *created_feature = created_features.Kth(j);
            R2Point created_feature_position(created_feature->scan_position.X(), created_feature->scan_position.Y());
            RNScalar distance_squared = R2SquaredDistance(local_maximum_position, created_feature_position);
            if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
          }

          // Create feature if not too close to previous (stronger) one
          if (well_spaced) {
            int grid_i = -1, grid_j = -1;
            dh_maxima_grid.IndexToIndices(grid_index, grid_i, grid_j);
            int scanline_index = (int) (dh_scanline_index_grid.GridValue(grid_index) + 0.5);
            int point_index = (int) (dh_point_index_grid.GridValue(grid_index) + 0.5);
            GSVScanline *scanline = scan->Scanline(scanline_index);
            RNScalar z = scanline->EstimatedGroundZ();
            if (z == RN_UNKNOWN) continue;
            R3Point position(x, y, z); 
            RNScalar score = dh_maxima_grid.GridValue(grid_index);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_POINT_FEATURE_TYPE, scanline, point_index, position, R3zero_vector, R3zero_vector, 1, score);
            ComputeSpinImageDescriptor(feature->descriptor, dh_position_x_grid, dh_position_y_grid, grid_i, grid_j);
            created_features.Insert(feature);
            InsertFeature(feature);
          }
        }

        // Release points
        if (!scan->ReleasePoints()) {
          fprintf(stderr, "Unable to release points for scan\n");
          return 0;
        }

#if 0
        // Write grids for debugging
        char buffer[1024];
        sprintf(buffer, "tmp/%d_%d_%d_dh_vertical_valley.grd", ir, is, ia);
        dh_vertical_valley_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_vote.grd", ir, is, ia);
        dh_vote_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_maxima.grd", ir, is, ia);
        dh_maxima_grid.WriteFile(buffer);
#endif
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Curb feature creation functions
////////////////////////////////////////////////////////////////////////

static int
ReadCurbFeatures(GSVPoseOptimization *optimization, GSVScan *scan, const char *filename)
{
  // Parameters
  const RNAngle max_curvature = 1;

  // Open curb path file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open curb path file %s\n", filename);
    return 0;
  }

  // Read XYZ positions of curb path vertices
  RNScalar score[3];
  R3Point position[3];
  int scanline_index[3], point_index[3], count;
  while (fscanf(fp, "%lf%lf%lf%lf%d%d%d\n", &position[2][0], &position[2][1], &position[2][2], 
    &score[2], &scanline_index[2], &point_index[2], &count) == (unsigned int) 7) {
    // Create feature
    if (count >= 2) {
      R3Vector va = position[0] - position[1];
      R3Vector vb = position[2] - position[1];
      RNLength lena = va.Length();
      RNLength lenb = vb.Length();
      if (RNIsNotZero(lena) && RNIsNotZero(lenb)) {
        RNScalar angle = R3InteriorAngle(va, vb);
        RNScalar curvature = (RN_PI - angle) / (lena + lenb);
        if (fabs(curvature) <= max_curvature) {
          R3Vector direction = position[2] - position[0];
          RNLength length = direction.Length();
          if (RNIsNotZero(length)) {
            RNScalar scale = 0.5 * length;
            direction /= length;
            R3Vector normal = direction % R3posz_vector;
            normal.Normalize();
            if (scan->SegmentIndex() == 0) normal.Flip();
            GSVScanline *scanline = scan->Scanline(scanline_index[1]);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index[1],  
              position[1], direction, normal, scale, score[1]);
            optimization->InsertFeature(feature);
          }
        }
      }
    }

    // Remember stuff
    position[0] = position[1];
    scanline_index[0] = scanline_index[1];
    point_index[0] = point_index[1];
    score[0] = score[1];
    position[1] = position[2];
    scanline_index[1] = scanline_index[2];
    point_index[1] = point_index[2];
    score[1] = score[2];
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateScanCurbFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/scan_curb_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Run program to extract curbs from scans
  char gsv2map_command[4096];
  sprintf(gsv2map_command, "gsv2map %s %s", scene->Filename(), tmp_directory);
  system(gsv2map_command);
    
  // Create curb features for every scan in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        char curb_name[1024];
        sprintf(curb_name, "%s/%02d_%02d_XYZ_CurbPath.txt", run_directory, is, ia);
        if (!ReadCurbFeatures(this, scan, curb_name)) return 0;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Edge feature creation functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
CreateScanEdgeFeatures(void)
{
  // Parameters
  const int min_gap_size = 5;
  const RNScalar min_height = 2;
  const RNScalar max_height = 20;
  const RNScalar min_depth_difference = 0.5 * RN_INFINITY; // boundaries with unknown
  const RNScalar expected_edge_height = 5.0;
  const RNScalar min_feature_spacing = 2;
  const RNScalar min_coverage = 0.2;
  const RNScalar depth_sigma = 1.0;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create edge features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;

        // Read points (needed for EstimatedGroundZ)
        if (!scan->ReadPoints()) {
          fprintf(stderr, "Unable to read points for scan\n");
          return 0;
        }

        // Read DH grids
        char image_name[4096];
        R2Grid dh_viewpoint_depth_grid;
        R2Grid dh_scanline_index_grid;
        R2Grid dh_point_index_grid;
        R2Grid dh_position_x_grid;
        R2Grid dh_position_y_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_ViewpointDepth.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_viewpoint_depth_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_y_grid.Read(image_name)) return 0;

        // Create DH grid containing 1 along vertical edges between large gaps of unknown and surface
        R2Grid dx_grid(dh_viewpoint_depth_grid);
        dx_grid.FillHoles((int) (dx_grid.WorldToGridScaleFactor() * min_gap_size));
        dx_grid.Substitute(R2_GRID_UNKNOWN_VALUE, RN_INFINITY);
        R2Grid dy_grid(dx_grid);
        dx_grid.Gradient(RN_X);
        dy_grid.Gradient(RN_Y);
        dx_grid.Abs();
        dy_grid.Abs();
        dx_grid.Threshold(min_depth_difference, 0, 1);
        dy_grid.Threshold(min_depth_difference, 0, 1);
        R2Grid dh_edge_grid(dx_grid);
        dh_edge_grid.Subtract(dy_grid);
        dh_edge_grid.Threshold(0.5, 0, 1);

        // Create DH vote grid
        R2Grid dh_vote_grid(dh_edge_grid);
        int y1 = (int) (dh_vote_grid.WorldToGridScaleFactor() * min_height);
        int y2 = (int) (dh_vote_grid.WorldToGridScaleFactor() * max_height);
        if (y1 >= dh_vote_grid.YResolution()) y1 = dh_vote_grid.YResolution()-1;
        if (y2 >= dh_vote_grid.YResolution()) y2 = dh_vote_grid.YResolution()-1;
        double y_sigma = 0.5 * dh_vote_grid.WorldToGridScaleFactor() * expected_edge_height;
        FilterVotes(dh_vote_grid, dh_viewpoint_depth_grid, y1, y2, y_sigma, depth_sigma);

        // Find local maxima of dh_vote_grid
        RNArray<const RNScalar *> local_maxima;
        R2Grid dh_maxima_grid(dh_vote_grid);
        dh_maxima_grid.MaskNonMaxima(0.25 * min_feature_spacing * dh_maxima_grid.WorldToGridScaleFactor());
        dh_maxima_grid.Threshold(min_coverage, R2_GRID_UNKNOWN_VALUE, R2_GRID_KEEP_VALUE);
        for (int i = 0; i < dh_maxima_grid.XResolution(); i++) {
          // Find maximum value in column
          int maximum_j = -1;
          RNScalar maximum_value = -FLT_MAX;
          for (int j = y1; j <= y2; j++) {
            RNScalar value = dh_maxima_grid.GridValue(i, j);
            if (value == R2_GRID_UNKNOWN_VALUE) continue;
            if (value == 0) continue;
            if (value < maximum_value) continue;
            maximum_value = value;
            maximum_j = j;
          }

          // Insert local maximum
          if (maximum_j >= 0) {
            const RNScalar *local_maximum = &dh_maxima_grid(i, maximum_j);
            local_maxima.Insert(local_maximum);
          }
        }

        // Sort local maxima
        local_maxima.Sort(CompareScalarPointers);

        // Create features at local maxima with minimum spacing
        RNArray<GSVFeature *> created_features;
        RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;
        for (int k = 0; k < local_maxima.NEntries(); k++) {
          // Get local maximum info
          const RNScalar *local_maximum = local_maxima.Kth(k);
          int grid_index = local_maximum - dh_maxima_grid.GridValues();
          assert((grid_index >= 0) && (grid_index < dh_maxima_grid.NEntries()));

          // Determine the x,y position
          RNScalar x = dh_position_x_grid.GridValue(grid_index);
          RNScalar y = dh_position_y_grid.GridValue(grid_index);

          // Check if there is a previously created (stronger) feature within minimum spacing
          RNBoolean well_spaced = TRUE;
          R2Point local_maximum_position(x, y);
          for (int j = 0; j < created_features.NEntries(); j++) {
            GSVFeature *created_feature = created_features.Kth(j);
            R2Point created_feature_position(created_feature->scan_position.X(), created_feature->scan_position.Y());
            RNScalar distance_squared = R2SquaredDistance(local_maximum_position, created_feature_position);
            if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
          }

          // Create feature if not too close to previous (stronger) one
          if (well_spaced) {
            int grid_i = -1, grid_j = -1;
            dh_maxima_grid.IndexToIndices(grid_index, grid_i, grid_j);
            int scanline_index = (int) (dh_scanline_index_grid.GridValue(grid_index) + 0.5);
            int point_index = (int) (dh_point_index_grid.GridValue(grid_index) + 0.5);
            GSVScanline *scanline = scan->Scanline(scanline_index);
            RNScalar z = scanline->EstimatedGroundZ();
            R3Point position(x, y, z);
            RNScalar score = dh_maxima_grid.GridValue(grid_index);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_POINT_FEATURE_TYPE, scanline, point_index, position, R3zero_vector, R3zero_vector, 1.0, score);
            ComputeSpinImageDescriptor(feature->descriptor, dh_position_x_grid, dh_position_y_grid, grid_i, grid_j);
            created_features.Insert(feature);
            InsertFeature(feature);
          }
        }

        // Release points
        if (!scan->ReleasePoints()) {
          fprintf(stderr, "Unable to release points for scan\n");
          return 0;
        }

#if 0
        // Write grids for debugging
        char buffer[1024];
        sprintf(buffer, "tmp/%d_%d_%d_dh_edge.grd", ir, is, ia);
        dh_edge_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_vote.grd", ir, is, ia);
        dh_vote_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_maxima.grd", ir, is, ia);
        dh_maxima_grid.WriteFile(buffer);
#endif
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Plane feature creation functions
////////////////////////////////////////////////////////////////////////

static RNScalar *
CreatePlaneSegmentAreas(int nsegments, const R2Grid& id_grid, 
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid)
{
  // Allocate areas
  RNScalar *areas = new RNScalar [ nsegments ];
  if (!areas) {
    fprintf(stderr, "Unable to allocate plane segment areas\n");
    return NULL;
  }

  // Initialize areas
  for (int i = 0; i < nsegments; i++) areas[i] = 0;

  // Compute areas
  for (int iy = 0; iy < id_grid.YResolution()-1; iy++) {
    for (int ix = 0; ix < id_grid.XResolution()-1; ix++) {
      // Get segment id
      RNScalar id0_value = id_grid.GridValue(ix, iy);
      if (id0_value == R2_GRID_UNKNOWN_VALUE) continue;
      int id0 = (int) (id0_value + 0.5);
      assert((id0 >= 0) && (id0 < nsegments));

      // Determine position
      RNScalar x0 = position_x_grid.GridValue(ix, iy);
      RNScalar y0 = position_y_grid.GridValue(ix, iy);
      RNScalar z0 = position_z_grid.GridValue(ix, iy);
      if (x0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (y0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (z0 == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point p0(x0, y0, z0);
      
      // Determine dx
      RNScalar dx = 0;
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix+1, iy);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;

        // Get position
        RNScalar x1 = position_x_grid.GridValue(ix+1, iy);
        RNScalar y1 = position_y_grid.GridValue(ix+1, iy);
        RNScalar z1 = position_z_grid.GridValue(ix+1, iy);
        if (x1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (y1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (z1 == R2_GRID_UNKNOWN_VALUE) continue;
        R3Point p1(x1, y1, z1);

        // Compute distance
        dx = R3Distance(p0, p1);
      }

      // Determine dy
      RNScalar dy = 0;
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix, iy+1);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;

        // Get position
        RNScalar x1 = position_x_grid.GridValue(ix, iy+1);
        RNScalar y1 = position_y_grid.GridValue(ix, iy+1);
        RNScalar z1 = position_z_grid.GridValue(ix, iy+1);
        if (x1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (y1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (z1 == R2_GRID_UNKNOWN_VALUE) continue;
        R3Point p1(x1, y1, z1);

        // Compute distance
        dy = R3Distance(p0, p1);
      }

      // Check diagonal point
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix+1, iy+1);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;
      }

      // Sum area of point sample
      areas[id0] += dx * dy;
    }
  }

  // Return areas
  return areas;
}



static RNScalar *
CreatePlaneSegmentConnectivities(int nsegments, const R2Grid& id_grid)
{
  // Allocate connectivies
  RNScalar *connectivities = new RNScalar [ nsegments ];
  if (!connectivities) {
    fprintf(stderr, "Unable to allocate plane segment connectivities\n");
    return NULL;
  }

  // Allocate temporary data
  int *counts = new int [ nsegments ];

  // Initialize counts
  for (int i = 0; i < nsegments; i++) {
    connectivities[i] = 0.0;
    counts[i] = 0;
  }

  // Compute counts of neighbors in same segment
  for (int iy = 0; iy < id_grid.YResolution()-1; iy++) {
    for (int ix = 0; ix < id_grid.XResolution()-1; ix++) {
      // Get segment id
      RNScalar id0_value = id_grid.GridValue(ix, iy);
      if (id0_value == R2_GRID_UNKNOWN_VALUE) continue;
      int id0 = (int) (id0_value + 0.5);
      assert((id0 >= 0) && (id0 < nsegments));
      for (int s = -1; s <= 1; s++) {
        for (int t = -1; t <= 1; t++) {
          if ((s == 0) && (t == 0)) continue;
          RNScalar id1_value = id_grid.GridValue(ix+s, iy+t);
          if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
          int id1 = (int) (id1_value + 0.5);
          assert((id1 >= 0) && (id1 < nsegments));
          if (id0 == id1) connectivities[id0] += 1.0;
        }
      }
      counts[id0] += 8;
    }
  }

  // Convert counts to fractions
  for (int i = 0; i < nsegments; i++) {
    if (counts[i] == 0) continue;
    connectivities[i] /= counts[i];
  }

  // Delete temporary data
  delete [] counts;
  
  // Return connectivities
  return connectivities;
}



int GSVPoseOptimization::
CreateScanPlaneFeatures(void)
{
  // Parameters
  const RNScalar min_feature_spacing = 5;
  const RNScalar min_segment_area = 10;
  const RNScalar min_segment_connectivity = 0.75;
  const RNScalar target_segment_area = 200;
  const int min_segment_nentries = 100;
  RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create plane features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;

        // Read DA grids
        char image_name[4096];
        R2Grid da_scanline_index_grid;
        R2Grid da_point_index_grid;
        R2Grid da_position_x_grid;
        R2Grid da_position_y_grid;
        R2Grid da_position_z_grid;
        R2Grid da_plane_segment_A_grid;
        R2Grid da_plane_segment_B_grid;
        R2Grid da_plane_segment_C_grid;
        R2Grid da_plane_segment_D_grid;
        R2Grid da_plane_segment_id_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_y_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentA.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_A_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentB.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_B_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentC.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_C_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentD.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_D_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PlaneSegmentId.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_id_grid.Read(image_name)) return 0;

        // Mask scan points on GSV car
        for (int iy = 0; iy < 25; iy++) {
          for (int ix = 0; ix < da_plane_segment_id_grid.XResolution(); ix++) {
            da_scanline_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_point_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_x_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_y_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_z_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_A_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_B_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_C_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_D_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_id_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
          }
        }

        // Allocate temporary segment data 
        int nsegments = (int) (da_plane_segment_id_grid.Maximum() + 0.5) + 1;
        RNScalar *segment_connectivities = CreatePlaneSegmentConnectivities(nsegments, da_plane_segment_id_grid); 
        RNScalar *segment_areas = CreatePlaneSegmentAreas(nsegments, da_plane_segment_id_grid, 
          da_position_x_grid, da_position_y_grid, da_position_z_grid);
        RNArray<const RNScalar *> *segment_entries = new RNArray<const RNScalar *> [ nsegments ];
        for (int i = 0; i < da_plane_segment_id_grid.NEntries(); i++) {
          RNScalar segment_id_value = da_plane_segment_id_grid.GridValue(i);
          if (segment_id_value == R2_GRID_UNKNOWN_VALUE) continue;
          int segment_id = (int) (segment_id_value + 0.5);
          assert((segment_id >= 0) && (segment_id < nsegments));
          const RNScalar *segment_entry = &da_plane_segment_id_grid(i);
          segment_entries[segment_id].Insert(segment_entry);
        }

        // Create features at points sampled from within plane segment
        for (int segment_id = 0; segment_id < nsegments; segment_id++) {
          // Check segment area
          RNScalar segment_area = segment_areas[segment_id];
          if (segment_area < min_segment_area) continue;

          // Check segment connectivity
          RNScalar segment_connectivity = segment_connectivities[segment_id];
          if (segment_connectivity < min_segment_connectivity) continue;
          
          // Check segment entries
          int segment_nentries = segment_entries[segment_id].NEntries();
          if (segment_nentries < min_segment_nentries) continue;

          // printf("%6d %6d %9.3g %9.3g\n", segment_id, segment_entries[segment_id].NEntries(), segment_area, segment_connectivity);

          // Create features on segment with minimum spacing
          RNArray<GSVFeature *> created_features;
          for (int i = 0; i < segment_entries[segment_id].NEntries(); i++) {
            int grid_ix, grid_iy;
            const RNScalar *segment_entry = segment_entries[segment_id].Kth(i);
            int grid_index = segment_entry - da_plane_segment_id_grid.GridValues();
            da_plane_segment_id_grid.IndexToIndices(grid_index, grid_ix, grid_iy);

            // Determine the position
            RNScalar x = da_position_x_grid.GridValue(grid_index);
            RNScalar y = da_position_y_grid.GridValue(grid_index);
            RNScalar z = da_position_z_grid.GridValue(grid_index);
            if (x == R2_GRID_UNKNOWN_VALUE) continue;
            if (y == R2_GRID_UNKNOWN_VALUE) continue;
            if (z == R2_GRID_UNKNOWN_VALUE) continue;
            R3Point position(x, y, z);

            // Check if there is a previously created feature within minimum spacing
            RNBoolean well_spaced = TRUE;
            for (int j = 0; j < created_features.NEntries(); j++) {
              GSVFeature *created_feature = created_features.Kth(j);
              RNScalar distance_squared = R3SquaredDistance(position, created_feature->scan_position);
              if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
            }

             // Create feature if not too close to previously created one
            if (well_spaced) {
              // Get feature information
              RNScalar scanline_index_value = da_scanline_index_grid.GridValue(grid_index);
              RNScalar point_index_value = da_point_index_grid.GridValue(grid_index);
              RNScalar a = da_plane_segment_A_grid.GridValue(grid_index);
              RNScalar b = da_plane_segment_B_grid.GridValue(grid_index);
              RNScalar c = da_plane_segment_C_grid.GridValue(grid_index);
              if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (a == R2_GRID_UNKNOWN_VALUE) continue;
              if (b == R2_GRID_UNKNOWN_VALUE) continue;
              if (c == R2_GRID_UNKNOWN_VALUE) continue;
              int scanline_index = (int) (scanline_index_value + 0.5);
              int point_index = (int) (point_index_value + 0.5);
              GSVScanline *scanline = scan->Scanline(scanline_index);
              R3Vector normal(a, b, c);
              normal.Normalize();
              RNScalar score = (segment_area < target_segment_area) ? segment_area / target_segment_area : 1.0;
              GSVFeature *feature = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index,  
                position, R3zero_vector, normal, min_feature_spacing, score);
              created_features.Insert(feature);
              InsertFeature(feature);
            }
          }

#if 0
          // Create correspondences between features within same plane segment
          if (created_features.NEntries() > 0) {
            double max_pairs = 30;
            double num_pairs = created_features.NEntries() * (created_features.NEntries() - 1) / 2;
            RNScalar pair_probability = max_pairs / num_pairs;
            for (int i = 0; i < created_features.NEntries(); i++) {
              GSVFeature *feature1 = created_features.Kth(i);
              if (feature1->scan_normal[2] > 0.5) continue; // hack to include only vertical surfaces
              for (int j = i+1; j < created_features.NEntries(); j++) {
                if (RNRandomScalar() > pair_probability) continue;
                GSVFeature *feature2 = created_features.Kth(j);
                GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature1, feature2, 0.99);
                InsertCorrespondence(correspondence);
              }
            }
          }
#endif
        }

        // Delete temporary segment data 
        delete [] segment_connectivities;
        delete [] segment_areas;
        delete [] segment_entries;
      }
    }
  }

  // Return success
  return 1;
}




////////////////////////////////////////////////////////////////////////
// PoseOptimization correspondence creation functions
////////////////////////////////////////////////////////////////////////

RNScalar GSVPoseOptimization::
CorrespondenceScore(const GSVFeature *feature0, const GSVFeature *feature1)
{
  // Check types
  if (feature0 == feature1) return -1;
  if (feature0->feature_type != feature1->feature_type) return -1;

  // Check scanlines
  GSVScanline *scanline0 = feature0->scanline;
  GSVScanline *scanline1 = feature1->scanline;
  if (scanline0 && scanline1 && (scanline0 == scanline1)) return -1;

  // Check images 
  GSVImage *image0 = feature0->image;
  GSVImage *image1 = feature1->image;
  if (image0 && image1 && (image0 == image1)) return -1;

  // Check path parameters
  if (min_intracorrespondence_path_parameter_difference > 0) {
    GSVScan *scan0 = scanline0->Scan();
    GSVScan *scan1 = scanline1->Scan();
    GSVSegment *segment0 = scan0->Segment();
    GSVSegment *segment1 = scan1->Segment();
    if (segment0 == segment1) {
      ScanlineData *scanline_data0 = (ScanlineData *) scanline0->Data();
      ScanlineData *scanline_data1 = (ScanlineData *) scanline1->Data();
      if (scanline_data0 && scanline_data1) {
        RNScalar u0 = scanline_data0->path_parameter;
        RNScalar u1 = scanline_data1->path_parameter;
        RNScalar path_parameter_difference = fabs(u0 - u1);
        if (path_parameter_difference < min_intracorrespondence_path_parameter_difference) return -1;
        if (path_parameter_difference < max_intracorrespondence_euclidean_distance) return -1; // ??? TEMPORARY ???
      }
    } 
  }

  // Check positions
  RNScalar euclidean_distance_score = 1;
  if (max_intracorrespondence_euclidean_distance > 0) {
    RNScalar max_intracorrespondence_euclidean_distance_squared = max_intracorrespondence_euclidean_distance * max_intracorrespondence_euclidean_distance;
    if (feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) max_intracorrespondence_euclidean_distance_squared *= 1.5;
    if (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) max_intracorrespondence_euclidean_distance_squared *= 1.5;
    R3Point position0 = FeaturePosition(feature0, TRUE);
    R3Point position1 = FeaturePosition(feature1, TRUE);
    RNScalar euclidean_distance_squared = R3SquaredDistance(position0, position1);
    if (euclidean_distance_squared > max_intracorrespondence_euclidean_distance_squared) return -1;
    euclidean_distance_score = exp(-euclidean_distance_squared / max_intracorrespondence_euclidean_distance_squared);
  }

  // Check normals
  RNScalar normal_angle_score = 1;
  if (max_intracorrespondence_normal_angle > 0) {
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      RNScalar max_intracorrespondence_normal_angle_squared = max_intracorrespondence_normal_angle * max_intracorrespondence_normal_angle;
      R3Vector normal0 = FeatureNormal(feature0, TRUE);
      R3Vector normal1 = FeatureNormal(feature1, TRUE);
      if (RNIsZero(normal0.Dot(normal0))) return -1;
      if (RNIsZero(normal1.Dot(normal1))) return -1;
      RNScalar dot = normal0.Dot(normal1);
      if (dot < 0) return -1;
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_intracorrespondence_normal_angle) return -1;
      RNScalar normal_angle_squared = normal_angle * normal_angle;
      normal_angle_score = exp(-normal_angle_squared / max_intracorrespondence_normal_angle_squared);
    }
  }

  // Check directions
  RNScalar direction_angle_score = 1;
  if (max_intracorrespondence_direction_angle > 0) {
    if ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) {
      RNScalar max_intracorrespondence_direction_angle_squared = max_intracorrespondence_direction_angle * max_intracorrespondence_direction_angle;
      R3Vector direction0 = FeatureDirection(feature0, TRUE);
      R3Vector direction1 = FeatureDirection(feature1, TRUE);
      if (RNIsZero(direction0.Dot(direction0))) return -1;
      if (RNIsZero(direction1.Dot(direction1))) return -1;
      RNScalar dot = direction0.Dot(direction1);
      if (dot < 0) return -1;
      RNAngle direction_angle = (dot < 1) ? acos(dot) : 0;
      if (direction_angle > max_intracorrespondence_direction_angle) return -1;
      RNScalar direction_angle_squared = direction_angle * direction_angle;
      direction_angle_score = exp(-direction_angle_squared / max_intracorrespondence_direction_angle_squared);
    }
  }

  // Check descriptors
  RNScalar descriptor_distance_score = 1;
  const GSVDescriptor& descriptor0 = feature0->descriptor;
  const GSVDescriptor& descriptor1 = feature1->descriptor;
  if (descriptor0.descriptor_type != descriptor1.descriptor_type) return -1;
  RNScalar max_d = 0;
  if (descriptor0.descriptor_type == GSV_SPIN_IMAGE_DESCRIPTOR_TYPE) 
    max_d = max_intracorrespondence_spin_image_descriptor_distance;
  else if (descriptor0.descriptor_type == GSV_SIFT_DESCRIPTOR_TYPE) 
    max_d = max_intracorrespondence_spin_image_descriptor_distance;
  else if (descriptor0.descriptor_type == GSV_LINE_DESCRIPTOR_TYPE) 
    max_d = max_intracorrespondence_line_descriptor_distance;
  RNScalar max_dd = max_d * max_d;
  if (max_dd > 0) {
    RNScalar dd = descriptor0.SquaredDistance(descriptor1);
    if (dd > max_dd) return -1;
    descriptor_distance_score = exp(-dd/max_dd);
  }

  // Compute score
  RNScalar score = feature0->score * feature1->score * euclidean_distance_score * 
    descriptor_distance_score * direction_angle_score * normal_angle_score;
  //  printf("%9.6f : %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", score, 
  //   feature0->score, feature1->score, euclidean_distance_score, 
  //   descriptor_distance_score, direction_angle_score, normal_angle_score);

  // Return score
  return score;
}



RNBoolean GSVPoseOptimization::
IsMatchFeasible(const GSVFeature *feature0, const GSVFeature *feature1, const GSVFeatureCorrespondenceSet *match)
{
  // Check if features are compatible
  RNScalar score = CorrespondenceScore(feature0, feature1);
  if (score <= 0) return FALSE;
  if (!match) return TRUE;

  // Useful variables
  GSVFeature *feature[2][2];
  GSVScanline *scanline[2][2];
  GSVScan *scan[2][2];
  GSVSegment *segment[2][2];
  R3Point position[2][2];
  feature[0][0] = (GSVFeature *) feature0;
  feature[0][1] = (GSVFeature *) feature1;
  scanline[0][0] = feature[0][0]->scanline;
  scanline[0][1] = feature[0][1]->scanline;
  scan[0][0] = scanline[0][0]->Scan();
  scan[0][1] = scanline[0][1]->Scan();
  segment[0][0] = scan[0][0]->Segment();
  segment[0][1] = scan[0][1]->Segment();
  position[0][0] = FeaturePosition(feature[0][0], TRUE);
  position[0][1] = FeaturePosition(feature[0][1], TRUE);

  // Check for compatibility with correspondences in match
  for (int i1 = 0; i1 < match->NCorrespondences(); i1++) {
    GSVFeatureCorrespondence *c1 = match->Correspondence(i1);
      
    // Get correspondence info 
    feature[1][0] = c1->Feature(0);
    feature[1][1] = c1->Feature(1);
    scanline[1][0] = feature[1][0]->scanline;
    scanline[1][1] = feature[1][1]->scanline;
    scan[1][0] = scanline[1][0]->Scan();
    scan[1][1] = scanline[1][1]->Scan();
    segment[1][0] = scan[1][0]->Segment();
    segment[1][1] = scan[1][1]->Segment();
    position[1][0] = FeaturePosition(feature[1][0], TRUE);
    position[1][1] = FeaturePosition(feature[1][1], TRUE);

    // Check pairs of pairs of features
    for (int j0 = 0; j0 < 2; j0++) {
      for (int j1 = 0; j1 < 2; j1++) {
        // Check if first pair of features are in same segment
        if (segment[0][j0] == segment[1][j1]) {
          // Check if second pair of features are in same segment
          if ((segment[0][1-j0] == segment[1][1-j1]) && (j0 == 0)) {
            // Check if correspondence is duplicate with another correspondence in match
            if ((feature[0][j0] == feature[1][j1]) && (feature[0][1-j0] == feature[1][1-j1])) return FALSE;

            // Check if feature is already corresponding to another feature in the same segment
            if ((feature[0][j0] == feature[1][j1]) && (feature[0][1-j0] != feature[1][1-j1])) return FALSE;

            // Check if path parameters of first pair of features are too close (?????)
            if (min_intrasegment_path_parameter_difference > 0) {
              ScanlineData *scanline_data0 = (ScanlineData *) scanline[0][j0];
              ScanlineData *scanline_data1 = (ScanlineData *) scanline[1][j1];
              if (scanline_data0 && scanline_data1) {
                RNScalar path_parameter_difference = fabs(scanline_data0->path_parameter - scanline_data1->path_parameter);
                if (path_parameter_difference < min_intrasegment_path_parameter_difference) return FALSE;
              }
            } 

            // Check ratio of distances between corresponding features in same segments
            if (max_intrasegment_distance_ratio > 0) {
              R3Point position10 = FeaturePosition(feature[1][j1], TRUE);
              R3Point position11 = FeaturePosition(feature[1][1-j1], TRUE);
              RNLength d0 = R3Distance(position[0][j0], position[1][j1]);
              RNLength d1 = R3Distance(position[0][1-j0], position[1][1-j1]);
              RNLength dmax = (d0 > d1) ? d0 : d1;
              RNScalar distance_ratio = (dmax > 0) ? fabs(d0 - d1) / dmax : 0;
              if (distance_ratio > max_intrasegment_distance_ratio) return FALSE;
            }
          }
        }
      }
    }
  }

  // Passed all tests
  return TRUE;
}



int GSVPoseOptimization::
CreateCoplanarCorrespondences(void)
{
  // Useful variables
  double max_offplane_distance = 1;
  double min_normal_angle_cosine = cos(RN_PI/8.0);
  double max_normal_z = cos(RN_PI/8.0);

  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature1 = Feature(i);
    if (feature1->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
    GSVScanline *scanline1 = feature1->scanline;
    if (!scanline1) continue;
    GSVScan *scan1 = scanline1->Scan();
    if (!scan1) continue;
    GSVSegment *segment1 = scan1->Segment();
    if (!segment1) continue;
    R3Point position1 = feature1->scan_position;
    const R3Vector& normal1 = feature1->scan_normal;
    if (RNIsZero(normal1.Dot(normal1))) continue;
    if (normal1.Z() > max_normal_z) continue;
    R3Plane plane1(position1, normal1);
    for (int j = i+1; j < NFeatures(); j++) {
      GSVFeature *feature2 = Feature(j);
      if (feature2->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
      GSVScanline *scanline2 = feature2->scanline;
      if (!scanline2) continue;
      GSVScan *scan2 = scanline2->Scan();
      if (!scan2) continue;
      GSVSegment *segment2 = scan2->Segment();
      if (!segment2) continue;
      if (segment1 != segment2) continue;
      R3Point position2 = feature2->scan_position;
      const R3Vector& normal2 = feature2->scan_normal;
      if (RNIsZero(normal2.Dot(normal2))) continue;
      if (normal2.Z() > max_normal_z) continue;
      R3Plane plane2(position2, normal2);
      if (fabs(normal1.Dot(normal2)) < min_normal_angle_cosine) continue;
      if (R3Distance(plane2, position1) > max_offplane_distance) continue;
      if (R3Distance(plane1, position2) > max_offplane_distance) continue;
      GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature1, feature2, 1);
      InsertCorrespondence(correspondence);
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageToScanCorrespondences(void)
{
  // Useful variables
  double max_distance = 1;

  RNArray<GSVFeatureCorrespondence *> tmp(correspondences.correspondences);
  for (int i = 0; i < tmp.NEntries(); i++) {
    GSVFeatureCorrespondence *correspondence = tmp.Kth(i);
    for (int j = 0; j < 2; j++) {
      GSVFeature *feature = correspondence->Feature(j);
      if (!feature) continue;
      GSVImage *image = feature->image;
      if (!image) continue;
      GSVTapestry *tapestry = image->Tapestry();
      if (!tapestry) continue;
      GSVSegment *segment = tapestry->Segment();
      if (!segment) continue;

      // Find closest mesh point
      GSVMesh *closest_mesh = NULL;
      RNLength closest_distance = FLT_MAX;
      R3Point closest_point = R3zero_point;
      R3Vector closest_normal = R3zero_vector;
      GSVMeshVertex *closest_vertex = NULL;
      for (int k = 0; k < segment->NScans(); k++) {
        GSVScan *scan = segment->Scan(k);
        if (scan->NScanlines() == 0) continue;
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;
        R3MeshIntersection closest;
        if (mesh->ClosestPoint(feature->scan_position, &closest, 0, max_distance)) {
          if ((!closest_mesh) || (closest.t < closest_distance)) {
            closest_mesh = mesh;
            closest_distance = closest.t;
            closest_point = closest.point;
            closest_normal = mesh->FaceNormal(closest.face);
            if (closest.type == R3_MESH_VERTEX_TYPE) closest_vertex = (GSVMeshVertex *) closest.vertex;
            else if (closest.type == R3_MESH_EDGE_TYPE) closest_vertex = (GSVMeshVertex *) mesh->VertexOnEdge(closest.edge, 0);
            else if (closest.type == R3_MESH_FACE_TYPE) closest_vertex = (GSVMeshVertex *) mesh->VertexOnFace(closest.face, 0);
            else RNAbort("Invalid closest point type\n");
          }
        }
      }   
      
      // Create correspondence
      if (closest_mesh) {
        GSVScanline *scanline = closest_mesh->VertexScanline(closest_vertex);
        int point_index = closest_mesh->VertexPointIndex(closest_vertex);
        GSVFeature *feature2 = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index, closest_point, R3zero_vector, closest_normal, 1, 1);
        InsertFeature(feature2);
        GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature, feature2, 1);
        InsertCorrespondence(correspondence);
      }
    }
  }
 
  // Return success
  return 1;
}



static R3Point 
FeaturePositionCallback(GSVFeature *feature, void *data)
{
  GSVPoseOptimization *optimization = (GSVPoseOptimization *) data;
  return optimization->FeaturePosition(feature, TRUE);
}



int GSVPoseOptimization::
CreateMutuallyClosestCorrespondences(void)
{
  // Create array of candidate features
  RNArray<GSVFeature *> features;
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    if (feature->scan_position == R3unknown_point) continue;
    features.Insert(feature);
  }

  // Create kdtree of features
  R3Kdtree<GSVFeature *> kdtree(features, FeaturePositionCallback, this);

  // Allocate temporary data
  double **closest_score = new double * [ NFeatures() ];
  GSVFeature ***closest_feature = new GSVFeature ** [ NFeatures() ];
  for (int i = 0; i < NFeatures(); i++) {
    closest_score[i] = new double [ segments.NEntries() ];
    closest_feature[i] = new GSVFeature * [ segments.NEntries() ];
    for (int j = 0; j < segments.NEntries(); j++) {
      closest_score[i][j] = -1;
      closest_feature[i][j] = NULL;
    }
  }

  // Find closest feature within every segment
  for (int i1 = 0; i1 < NFeatures(); i1++) {
    GSVFeature *feature1 = Feature(i1);
    GSVScanline *scanline1 = feature1->scanline;
    R3Point position1 = FeaturePosition(feature1, TRUE);
    RNArray<GSVFeature *> nearby_features;
    if (!kdtree.FindAll(position1, 0, max_intracorrespondence_euclidean_distance, nearby_features)) continue;
    for (int k = 0; k < nearby_features.NEntries(); k++) {
      GSVFeature *feature2 = nearby_features.Kth(k);
      if (feature2 == feature1) continue;
      GSVScanline *scanline2 = feature2->scanline;
      if (scanline2 == scanline1) continue;
      GSVScan *scan2 = scanline2->Scan();
      if (!scan2) continue;
      GSVSegment *segment2 = scan2->Segment();
      if (!segment2) continue;
      SegmentData *segment_data2 = (SegmentData *) segment2->Data();
      if (!segment_data2) continue;
      int segment_index2 = segment_data2->index;
      RNLength score = CorrespondenceScore(feature1, feature2);
      if (score < 0) continue;
      if (score > closest_score[i1][segment_index2]) {
        closest_feature[i1][segment_index2] = feature2;
        closest_score[i1][segment_index2] = score;
      }
    }
  }

  // Create correspondences between mutually closest features
  for (int i1 = 0; i1 < NFeatures(); i1++) {
    GSVFeature *feature1 = Feature(i1);
    GSVScanline *scanline1 = feature1->scanline;
    GSVScan *scan1 = scanline1->Scan();
    GSVSegment *segment1 = scan1->Segment();
    int j1 = segment1->SceneIndex();
    for (int j2 = 0; j2 < segments.NEntries(); j2++) {
      // GSVSegment *segment2 = segments.Kth(j2);
      // if (segment1 == segment2) continue;
      RNScalar score = closest_score[i1][j2];
      if (score <= 0) continue;
      GSVFeature *feature2 = closest_feature[i1][j2];
      if (!feature2) continue;
      int i2 = feature2->index;
      if (i2 < 0) continue;
      if (closest_feature[i2][j1] != feature1) continue;
      GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature1, feature2, score);
      InsertCorrespondence(correspondence);
    }
  }
  
  // Delete temporary data
  for (int i = 0; i < NFeatures(); i++) {
    delete [] closest_score[i];
    delete [] closest_feature[i];
  }
  delete [] closest_score;
  delete [] closest_feature;
  
  // Return success
  return 1;
}
  


int GSVPoseOptimization::
CreateAllFeasibleCorrespondences(void)
{
  // Consider every pair of features
  for (int i0 = 0; i0 < NFeatures(); i0++) {
    GSVFeature *feature0 = Feature(i0);
    for (int i1 = 0; i1 < NFeatures(); i1++) {
      GSVFeature *feature1 = Feature(i1);
      RNScalar score = CorrespondenceScore(feature0, feature1);
      if (score > 0) {
        GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, score);
        InsertCorrespondence(correspondence);
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExhaustiveSearch(const GSVFeatureCorrespondenceSet& pairwise_matches, 
  int max_correspondences, GSVFeatureCorrespondenceSet& best_match)
{
  // Compute alignment and score
  Solve(TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE);
  correspondences.score = Score();
  if (correspondences.score > best_match.score) {
    printf("%d %d : %g\n", NFeatures(), NCorrespondences(), correspondences.score);
    best_match = correspondences;
  }

  // Check current number of correspondences
  if (correspondences.NCorrespondences() < max_correspondences) {
    for (int i = 0; i < pairwise_matches.NCorrespondences(); i++) {
      GSVFeatureCorrespondence *correspondence = pairwise_matches.Correspondence(i);
      InsertCorrespondence(correspondence);
      ExhaustiveSearch(pairwise_matches, max_correspondences, best_match);
      RemoveCorrespondence(correspondence);
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExhaustiveSearch(void)
{
  // Parameters
  int max_correspondences = 3;
  int num_initial_correspondences = NCorrespondences();

  // Create pairwise corrspondences
  if (!CreateAllFeasibleCorrespondences()) return 0;
  GSVFeatureCorrespondenceSet pairwise_matches(correspondences);
  correspondences.correspondences.Truncate(num_initial_correspondences);

  // Find best set of correspondences with branch-and-bound
  GSVFeatureCorrespondenceSet best_match(scene);
  if (!ExhaustiveSearch(pairwise_matches, max_correspondences, best_match)) {
    fprintf(stderr, "Branch and bound failed\n");
    return 0;
  }

  // Copy best match into correspondences
  correspondences = best_match;

  // Solve fully
  Solve();

  // Return success
  return 1;
}



int GSVPoseOptimization::
RANSAC(void)
{
  // Parameters
  int max_correspondences = 5;
  int max_iterations = 100;
  int num_initial_correspondences = NCorrespondences();

  // Create pairwise corrspondences
  if (!CreateAllFeasibleCorrespondences()) return 0;
  GSVFeatureCorrespondenceSet pairwise_matches(correspondences);
  correspondences.correspondences.Truncate(num_initial_correspondences);
  pairwise_matches.Sort();

  // Sample sets of correspondences 
  GSVFeatureCorrespondenceSet best_match(scene);
  for (int iter = 0; iter < max_iterations; iter++) {
    correspondences.correspondences.Truncate(num_initial_correspondences);
    for (int i = 0; i < max_correspondences; i++) {
      int index = (int) (RNRandomScalar() * pairwise_matches.NCorrespondences());
      GSVFeatureCorrespondence *correspondence = pairwise_matches.Correspondence(index);
      GSVFeature *feature0 = correspondence->Feature(0);
      GSVFeature *feature1 = correspondence->Feature(1);
      if (!IsMatchFeasible(feature0, feature1, &correspondences)) continue;
      InsertCorrespondence(correspondence);
    }
    if (NCorrespondences() > 0) {
      Solve(TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, 1.0);
      correspondences.score = Score();
      printf("%g\n", correspondences.score);
      if (correspondences.score > best_match.score) {
        best_match = correspondences;
      }
    }
  }

  // Copy best match into correspondences
  correspondences = best_match;

  // Solve fully
  Solve();

  // Return success
  return 1;
}



int GSVPoseOptimization::
ICP(void)
{
  // Just checking
  if (icp_max_distance_end <= 0) icp_max_distance_end = icp_max_distance_start;
  if (icp_max_distance_end >= icp_max_distance_start) icp_max_distance_end = icp_max_distance_start;

  // Remember parameters
  RNScalar saved = max_intracorrespondence_euclidean_distance;
  int num_initial_correspondences = NCorrespondences();

  // Iteratively create mutually closest correspondences
  RNScalar rigidity = 1;
  const int max_iterations = 32;
  const RNScalar min_rigidity = 0.25;
  max_intracorrespondence_euclidean_distance = icp_max_distance_start;
  for (int i = 0; i < max_iterations; i++) {
    correspondences.correspondences.Truncate(num_initial_correspondences);
    CreateMutuallyClosestCorrespondences();
    Solve(TRUE, FALSE,  FALSE, FALSE, TRUE, TRUE,  TRUE, TRUE, TRUE, rigidity);
    RNScalar score = Score();
    printf("%d : %d %d : %g : %g\n", i+1, NFeatures(), NCorrespondences(), 
      max_intracorrespondence_euclidean_distance, score);
    if (max_intracorrespondence_euclidean_distance == icp_max_distance_end) break;
    rigidity *= 0.5;
    if (rigidity < min_rigidity) {
      rigidity = min_rigidity;
    }
    max_intracorrespondence_euclidean_distance *= 0.5;
    if (max_intracorrespondence_euclidean_distance < icp_max_distance_end) {
      max_intracorrespondence_euclidean_distance = icp_max_distance_end;
    }
  }    

  // Solve fully
  Solve();

  // Restore original parameters
  max_intracorrespondence_euclidean_distance = saved;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateCorrespondences(void)
{
  // Use exhaustive search
  // ExhaustiveSearch();

  // Use RANSAC
  // RANSAC();

  // Use ICP
  ICP();

  // Return success
  return 1;
}



#if 0

RNScalar GSVPoseOptimization::
ScoreCorrespondence(GSVFeature *featureA, GSVFeature *featureB)
{
  // Parameters
  RNScalar similarity_sigma = 1;
  RNScalar similarity_scale = -1.0 / (2 * similarity_sigma * similarity_sigma);
  RNScalar max_distance = 3 * similarity_sigma;
  RNScalar max_squared_distance = max_distance * max_distance;

  // Get translation vector
  R3Point positionA = FeaturePosition(featureA, FALSE);
  R3Point positionB = FeaturePosition(featureB, FALSE);
  R3Vector translation = positionB - positionA;

  // Get segments
  GSVScanline *scanlineA = featureA->scanline;
  GSVScanline *scanlineB = featureB->scanline;
  GSVScanline *scanA = scanlineA->Scan();
  GSVScanline *scanB = scanlineB->Scan();
  GSVScanline *segmentA = scanA->Segment();
  GSVScanline *segmentB = scanB->Segment();

  // Initialize result
  RNScalar similarity = 0;
  int feature_count = 0;

  // Consider every feature in segmentA
  for (int i1 = 0; i1 < NFeatures(); i1++) {
    GSVFeature *feature1 = Feature(i1);
    if (feature1->feature_type != GSV_SCAN_POINT_FEATURE_TYPE) continue;
    GSVScanline *scanline1 = feature1->scanline;
    GSVScan *scan1 = scanline1->Scan();
    GSVSegment *segment1 = scan1->Segment();
    if (segment1 != segmentA) continue;
    R3Point position1 = FeaturePosition(feature1, TRUE);
    position1 += translation;
    feature_count++;

    // Find closest feature in segmentB
    RNLength closest_squared_distance = max_squared_distance;
    for (int i2 = 0; i2 < NFeatures(); i2++) {
      GSVFeature *feature2 = Feature(i2);
      if (feature2->feature_type != GSV_SCAN_POINT_FEATURE_TYPE) continue;
      GSVScanline *scanline2 = feature2->scanline;
      GSVScan *scan2 = scanline2->Scan();
      GSVSegment *segment2 = scan2->Segment();
      if (segment2 != segmentB) continue;
      R3Point position2 = FeaturePosition(feature2, TRUE);
      RNLength squared_distance = R3SquaredDistance(position1, position2);
      if (squared_distance < closest_squared_distance) {
        closest_squared_distance = squared_distance;
      }
    }

    // Add similarity
    if (closest_squared_distance < max_squared_distance) {
      similarity += exp(similarity_scale * closest_squared_distance);
    }
  }

  // Score is total similarity normalized by number of features
  RNScalar score = (feature_count > 0) ? similarity / feature_count : 0;
      
  // Return score
  return score;
}

#endif



RNScalar GSVPoseOptimization::
Score(void)
{
  // Parameters
  RNScalar similarity_sigma = 1;

  // Gather statistics
  RNScalar similarity = 0;
  RNScalar max_squared_distance = 0;
  RNScalar similarity_scale = -1.0 / (2 * similarity_sigma * similarity_sigma);;
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    assert(feature0 && feature1);
    GSVScanline *scanline0 = feature0->scanline;
    GSVScanline *scanline1 = feature1->scanline;
    assert(scanline0 && scanline1);
    R3Affine transformation0 = OptimizedTransformation(scanline0);
    R3Affine transformation1 = OptimizedTransformation(scanline1);
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      R3Point point0 = feature0->scan_position;
      R3Point point1 = feature1->scan_position;
      point0.Transform(transformation0);
      point1.Transform(transformation1);
      RNScalar squared_distance = R3SquaredDistance(point0, point1);
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }
    else if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      R3Plane plane0(feature0->scan_position, feature0->scan_normal);
      R3Point point1 = feature1->scan_position;
      plane0.Transform(transformation0);
      point1.Transform(transformation1);
      RNScalar distance = R3Distance(plane0, point1);
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
     }      
    else if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      R3Point point0 = feature0->scan_position;
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      point0.Transform(transformation0);
      plane1.Transform(transformation1);
      RNScalar distance = R3Distance(plane1, point0);
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }      
    else if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      R3Point point0 = feature0->scan_position;
      R3Point point1 = feature1->scan_position;
      R3Plane plane0(feature0->scan_position, feature0->scan_normal);
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      point0.Transform(transformation0);
      point1.Transform(transformation1);
      plane0.Transform(transformation0);
      plane1.Transform(transformation1);
      RNScalar distance0 = R3Distance(plane0, point1);
      RNScalar distance1 = R3Distance(plane1, point0);
      RNScalar distance = (distance0 > distance1) ? distance0 : distance1;
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }      
    else if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      R3Point point0 = feature0->scan_position;
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      point0.Transform(transformation0);
      ray1.Transform(transformation1);
      RNLength distance = R3Distance(ray1, point0);
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }      
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Point point1 = feature1->scan_position;
      ray0.Transform(transformation0);
      point1.Transform(transformation1);
      RNLength distance = R3Distance(ray0, point1);
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }      
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      ray0.Transform(transformation0);
      ray1.Transform(transformation1);
      // RNLength distance = R3Distance(ray0, ray1);
      RNLength distance = R3Distance(ray0.Line(), ray1.Line());
      RNScalar squared_distance = distance * distance;
      if (squared_distance > max_squared_distance) max_squared_distance = squared_distance;
      similarity += exp(similarity_scale * squared_distance);
    }      
  }

  // Score is total similarity normalized by number of features
  RNScalar score = (NFeatures() > 0) ? similarity / NFeatures() : 0;
      
  // Return score
  return score;
}



RNScalar GSVPoseOptimization::
RMSD(void)
{
  // Compute RMSD
  int count = 0;
  RNScalar sum = 0;
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    assert(feature0 && feature1);
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVScanline *scanline1 = feature1->scanline;
      assert(scanline0 && scanline1);
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Point point0 = feature0->scan_position;
      R3Point point1 = feature1->scan_position;
      point0.Transform(transformation0);
      point1.Transform(transformation1);
      sum += R3SquaredDistance(point0, point1);
      count++;
    }
    else if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVScanline *scanline1 = feature1->scanline;
      assert(scanline0 && scanline1);
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Plane plane0(feature0->scan_position, feature0->scan_normal);
      R3Point point1 = feature1->scan_position;
      plane0.Transform(transformation0);
      point1.Transform(transformation1);
      RNScalar distance = R3Distance(plane0, point1);
      sum += distance * distance;
      count++;
    }      
    else if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVScanline *scanline1 = feature1->scanline;
      assert(scanline0 && scanline1);
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Point point0 = feature0->scan_position;
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      point0.Transform(transformation0);
      plane1.Transform(transformation1);
      RNScalar distance = R3Distance(plane1, point0);
      sum += distance * distance;
      count++;
    }      
    else if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVScanline *scanline1 = feature1->scanline;
      assert(scanline0 && scanline1);
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Point point0 = feature0->scan_position;
      R3Point point1 = feature1->scan_position;
      R3Plane plane0(feature0->scan_position, feature0->scan_normal);
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      point0.Transform(transformation0);
      point1.Transform(transformation1);
      plane0.Transform(transformation0);
      plane1.Transform(transformation1);
      RNScalar distance0 = R3Distance(plane0, point1);
      RNScalar distance1 = R3Distance(plane1, point0);
      RNScalar distance = (distance0 > distance1) ? distance0 : distance1;
      sum += distance * distance;
      count++;
    }
    else if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVImage *image1 = feature1->image;
      assert(scanline0 && image1);
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(image1);
      R3Point point0 = feature0->scan_position;
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      point0.Transform(transformation0);
      ray1.Transform(transformation1);
      RNLength distance = R3Distance(ray1, point0);
      sum += distance * distance;
      count++;
    }      
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      GSVScanline *scanline1 = feature1->scanline;
      assert(image0 && scanline1);
      R3Affine transformation0 = OptimizedTransformation(image0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Point point1 = feature1->scan_position;
      ray0.Transform(transformation0);
      point1.Transform(transformation1);
      RNLength distance = R3Distance(ray0, point1);
      sum += distance * distance;
      count++;
    }      
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      GSVImage *image1 = feature1->image;
      assert(image0 && image1);
      R3Affine transformation0 = OptimizedTransformation(image0);
      R3Affine transformation1 = OptimizedTransformation(image1);
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      ray0.Transform(transformation0);
      ray1.Transform(transformation1);
      // RNLength distance = R3Distance(ray0, ray1);
      RNLength distance = R3Distance(ray0.Line(), ray1.Line());
      sum += distance * distance;
      count++;
    }      
  }

  // Return RMSD
  return (count > 0) ? sqrt(sum/count) : 0;
}



////////////////////////////////////////////////////////////////////////
// Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".txt", 4)) {
    if (!ReadAsciiFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".fet", 4)) {
    if (!ReadBinaryFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".sft", 4)) {
    if (!ReadSiftFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".lin", 4)) {
    if (!ReadLineFile(filename)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".txt", 4)) {
    if (!WriteAsciiFile(filename, apply_pose_transformations)) return 0;
  }
  else if (!strncmp(extension, ".fet", 4)) {
    if (!WriteBinaryFile(filename, apply_pose_transformations)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Ascii Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadAsciiFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadAscii(fp)) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteAsciiFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteAscii(fp, apply_pose_transformations)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
ReadAscii(FILE *fp)
{
  // Read file
  char keyword[1024];
  RNArray<GSVFeature *> read_features;
  while (fscanf(fp, "%s", keyword) == (unsigned int) 1) {
    if (!strcmp(keyword, "F")) {
      // Parse feature type
      int feature_type;
      if (fscanf(fp, "%d", &feature_type) != (unsigned int) 1) {
        fprintf(stderr, "Invalid feature command\n");
        return 0;
      }

      // Parse feature info
      if ((feature_type == GSV_SCAN_POINT_FEATURE_TYPE) ||
          (feature_type == GSV_SCAN_LINE_FEATURE_TYPE) ||
          (feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
        // Read scan info
        char run_name[1024];
        int segment_index, scan_index, scanline_index, point_index;
        double score, px, py, pz, dx, dy, dz;
        if (fscanf(fp, "%s%d%d%d%d%lf%lf%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &scan_index, &scanline_index, &point_index, &score,
          &px, &py, &pz, &dx, &dy, &dz) != (unsigned int) 12) {
          fprintf(stderr, "Error parsing feature %d\n", NFeatures());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (scan_index >= 0) {
          if (scan_index >= segment->NScans()) { 
            fprintf(stderr, "Bad scan %d\n", scan_index);
            return 0;
          }
        }
        GSVScan *scan = segment->Scan(scan_index);
        if (scanline_index >= scan->NScanlines()) {
          fprintf(stderr, "Bad scanline %d\n", scanline_index);
          return 0;
        }
        GSVScanline *scanline = scan->Scanline(scanline_index);
     
        // Create feature
        R3Point position(px, py, pz);
        R3Vector direction(dx, dy, dz);
        GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, direction, direction, 1.0, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else if ((feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
               (feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) {
        // Read image info
        char run_name[1024];
        int segment_index, panorama_index, image_index;
        double score, px, py, dx, dy;
        if (fscanf(fp, "%s%d%d%d%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &panorama_index, &image_index, &score,
          &px, &py, &dx, &dy) != (unsigned int) 9) {
          fprintf(stderr, "Error parsing feature %d\n", read_features.NEntries());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (panorama_index >= 0) {
          if (panorama_index >= segment->NPanoramas()) { 
            fprintf(stderr, "Bad panorama %d\n", panorama_index);
            return 0;
          }
        }
        GSVPanorama *panorama = segment->Panorama(panorama_index);
        if (image_index >= panorama->NImages()) {
          fprintf(stderr, "Bad image %d\n", image_index);
          return 0;
        }
        GSVImage *image = panorama->Image(image_index);

        // Create feature
        R2Point position(px, py);
        R2Vector direction(dx, dy);
        GSVFeature *feature = new GSVFeature(feature_type, image, position, direction, 1.0, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else {
        fprintf(stderr, "Unrecognized feature type %d\n", feature_type);
        return 0;
      }
    }
    else if (!strcmp(keyword, "D")) {
      // Read descriptor info
      int feature_index, descriptor_type, nvalues;
      if (fscanf(fp, "%d%d%d", &feature_index, &descriptor_type, &nvalues) != (unsigned int) 3) {
        fprintf(stderr, "Unable to parse descriptor\n");
        return 0;
      }

      // Check feature index
      if ((feature_index < 0) || (feature_index >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index);
        return 0;
      }

      // Read descriptor
      GSVFeature *feature = read_features.Kth(feature_index);
      GSVDescriptor& descriptor = feature->descriptor;
      descriptor.descriptor_type = descriptor_type;
      descriptor.nvalues = nvalues;
      descriptor.values = new RNScalar [ nvalues ];
      for (int i = 0; i < nvalues; i++) {
        if (fscanf(fp, "%lf", &descriptor.values[i]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to parse descriptor\n");
          return 0;
        }
      }
    }
    else if (!strcmp(keyword, "C")) {
      // Read correspondence
      int feature_index0, feature_index1;
      double score;
      if (fscanf(fp, "%d%d%lf", &feature_index0, &feature_index1, &score) != (unsigned int) 3) {
        fprintf(stderr, "Unable to parse correspondence %d\n", NCorrespondences());
        return 0;
      }


      // Check feature indices
      if ((feature_index0 < 0) || (feature_index0 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index0);
        return 0;
      }
      if ((feature_index1 < 0) || (feature_index1 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index1);
        return 0;
      }

      // Create correspondence
      GSVFeature *feature0 = read_features.Kth(feature_index0);
      GSVFeature *feature1 = read_features.Kth(feature_index1);
      GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, score);
      InsertCorrespondence(correspondence);
    }
    else if (!strcmp(keyword, "include")) {
      // Read filename
      char filename[2048];
      if (fscanf(fp, "%s", filename) != (unsigned int) 1) {
        fprintf(stderr, "Unable to parse include filename\n");
        return 0;
      }

      // Read included file
      if (!ReadFile(filename)) return 0;
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteAscii(FILE *fp, RNBoolean apply_pose_transformations) const
{
  // Write features
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);

    // Write feature type
    fprintf(fp, "F %2d  ", feature->feature_type);

    // Write feature info
    if ((feature->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) ||
        (feature->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) ||
        (feature->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Get scene info
      GSVScanline *scanline = feature->scanline;
      GSVScan *scan = (scanline) ? scanline->Scan() : NULL;
      GSVSegment *segment = (scan) ? scan->Segment() : NULL;
      GSVRun *run = (segment) ? segment->Run() : NULL;
      int scanline_index = (scanline) ? scanline->ScanIndex() : -1;
      int scan_index = (scan) ? scan->SegmentIndex() : -1;
      int segment_index = (segment) ? segment->RunIndex() : -1;
      const char *run_name = (run) ? run->Name() : "None";

      // Compute scan position
      R3Point scan_position = FeaturePosition(feature, apply_pose_transformations); 
      R3Vector scan_normal = FeatureNormal(feature, apply_pose_transformations); 

      // Write feature info
      fprintf(fp, "%30s %2d %2d %6d %3d   %8.6f   %12.6f %12.6f %12.6f  %8.6f %8.6f %8.6f\n", 
        run_name, segment_index, scan_index, scanline_index, feature->scan_point_index, feature->score,
        scan_position.X(), scan_position.Y(), scan_position.Z(), 
        scan_normal.X(), scan_normal.Y(), scan_normal.Z());
    }
    else if ((feature->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
             (feature->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) {
      // Get scene info
      GSVImage *image = feature->image;
      GSVPanorama *panorama = (image) ? image->Panorama() : NULL;
      GSVSegment *segment = (panorama) ? panorama->Segment() : NULL;
      GSVRun *run = (segment) ? segment->Run() : NULL;
      int image_index = (image) ? image->PanoramaIndex() : -1;
      int panorama_index = (panorama) ? panorama->SegmentIndex() : -1;
      int segment_index = (segment) ? segment->RunIndex() : -1;
      const char *run_name = (run) ? run->Name() : "None";

      // Write feature info
      fprintf(fp, "%30s %2d %6d %2d   %8.6f   %12.6f %12.6f   %8.6f %8.6f\n", 
        run_name, segment_index,  panorama_index, image_index, feature->score,
        feature->image_position.X(), feature->image_position.Y(), 
        feature->image_direction.X(), feature->image_direction.Y());
    }
  }

  ////////////////////////////////////////
  
  // Write descriptors
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    GSVDescriptor& descriptor = feature->descriptor;
    fprintf(fp, "D %6d %2d %3d  ", i, descriptor.descriptor_type, descriptor.nvalues);
    for (int i = 0; i < descriptor.nvalues; i++) {
      fprintf(fp, "%12.6f ", descriptor.values[i]);
    }
    fprintf(fp, "\n");
  }

  ////////////////////////////////////////
  
  // Write correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    fprintf(fp, "C %d %d %g\n", feature0->index, feature1->index, correspondence->score);
  }
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Binary Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadBinaryFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadBinary(fp)) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteBinaryFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteBinary(fp, apply_pose_transformations)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
ReadBinary(FILE *fp)
{
  // Read endian test
  int endian_test;
  if (fread(&endian_test, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to read correspondence file\n");
    return 0;
  }
  
  // Check endian test
  if (endian_test != 1) {
    fprintf(stderr, "Unable to read format of correspondence file\n");
    return 0;
  }

  // Read rest of header
  int dummy = 0;
  int major_version, minor_version;
  int nfeatures, ncorrespondences;
  fread(&major_version, sizeof(int), 1, fp);
  fread(&minor_version, sizeof(int), 1, fp);
  fread(&nfeatures, sizeof(int), 1, fp);
  fread(&ncorrespondences, sizeof(int), 1, fp);
  for (int i = 0; i < 59; i++) fread(&dummy, sizeof(int), 1, fp);

  // Read features
  for (int i = 0; i < nfeatures; i++) {
    GSVFeature *feature = new GSVFeature();
    assert(feature);
    int scanline_index, image_index, descriptor_nvalues;
    fread(&feature->feature_type, sizeof(int), 1, fp);
    fread(&image_index, sizeof(int), 1, fp);
    fread(&scanline_index, sizeof(int), 1, fp);
    fread(&feature->scan_point_index, sizeof(int), 1, fp);
    fread(&feature->scan_position, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_direction, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_normal, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_scale, sizeof(RNScalar), 1, fp);
    fread(&feature->image_position, sizeof(RNCoord), 2, fp);
    fread(&feature->image_direction, sizeof(RNCoord), 2, fp);
    fread(&feature->image_scale, sizeof(RNScalar), 1, fp);
    fread(&feature->image_t, sizeof(RNScalar), 1, fp);
    fread(&feature->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 21; j++) fread(&dummy, sizeof(int), 1, fp);
    fread(&feature->descriptor.descriptor_type, sizeof(int), 1, fp);
    fread(&descriptor_nvalues, sizeof(int), 1, fp);
    if (descriptor_nvalues > 0) {
      feature->descriptor.nvalues = descriptor_nvalues;
      feature->descriptor.values = new RNScalar [ feature->descriptor.nvalues ];
      fread(feature->descriptor.values, sizeof(RNScalar), feature->descriptor.nvalues, fp);
    }
    feature->scanline = (scanline_index >= 0) ? scanlines.Kth(scanline_index) : NULL;
    feature->image = (image_index >= 0) ? images.Kth(image_index) : NULL;
    InsertFeature(feature);
  }

  // Read correspondences
  for (int i = 0; i < ncorrespondences; i++) {
    GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence();
    assert(correspondence);
    int feature_index0, feature_index1;
    fread(&feature_index0, sizeof(int), 1, fp);
    fread(&feature_index1, sizeof(int), 1, fp);
    fread(&correspondence->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 13; j++) fread(&dummy, sizeof(int), 1, fp);
    correspondence->features[0] = (feature_index0 >= 0) ? Feature(feature_index0) : NULL;
    correspondence->features[1] = (feature_index1 >= 0) ? Feature(feature_index1) : NULL;
    InsertCorrespondence(correspondence);
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteBinary(FILE *fp, RNBoolean apply_pose_transformations) const
{
  // Convenient variables
  int dummy = 0;

  // Write endian test
  int endian_test = 1;
  if (fwrite(&endian_test, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to write correspondence file\n");
    return 0;
  }

  // Write rest of header
  int major_version = 0;
  int minor_version = 1;
  int nfeatures = NFeatures();
  int ncorrespondences = NCorrespondences();
  fwrite(&major_version, sizeof(int), 1, fp);
  fwrite(&minor_version, sizeof(int), 1, fp);
  fwrite(&nfeatures, sizeof(int), 1, fp);
  fwrite(&ncorrespondences, sizeof(int), 1, fp);
  for (int i = 0; i < 59; i++) fwrite(&dummy, sizeof(int), 1, fp);

  // Write features
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    GSVScanline *scanline = feature->scanline;
    ScanlineData *scanline_data = (scanline) ? (ScanlineData *) scanline->Data() : NULL;
    int scanline_index = (scanline_data) ? scanline_data->index : -1;
    R3Point scan_position = FeaturePosition(feature, apply_pose_transformations); 
    R3Vector scan_direction = FeatureDirection(feature, apply_pose_transformations); 
    R3Vector scan_normal = FeatureNormal(feature, apply_pose_transformations); 
    GSVImage *image = feature->image;
    ImageData *image_data = (image) ? (ImageData *) image->Data() : NULL;
    int image_index = (image_data) ? image_data->index : -1;
    const GSVDescriptor& descriptor = feature->descriptor;
    fwrite(&feature->feature_type, sizeof(int), 1, fp);
    fwrite(&image_index, sizeof(int), 1, fp);
    fwrite(&scanline_index, sizeof(int), 1, fp);
    fwrite(&feature->scan_point_index, sizeof(int), 1, fp);
    fwrite(&scan_position, sizeof(RNCoord), 3, fp);
    fwrite(&scan_direction, sizeof(RNCoord), 3, fp);
    fwrite(&scan_normal, sizeof(RNCoord), 3, fp);
    fwrite(&feature->scan_scale, sizeof(RNScalar), 1, fp);
    fwrite(&feature->image_position, sizeof(RNCoord), 2, fp);
    fwrite(&feature->image_direction, sizeof(RNCoord), 2, fp);
    fwrite(&feature->image_scale, sizeof(RNScalar), 1, fp);
    fwrite(&feature->image_t, sizeof(RNScalar), 1, fp);
    fwrite(&feature->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 21; j++) fwrite(&dummy, sizeof(int), 1, fp);
    fwrite(&descriptor.descriptor_type, sizeof(int), 1, fp);
    fwrite(&descriptor.nvalues, sizeof(int), 1, fp);
    if (descriptor.nvalues > 0) {
      fwrite(descriptor.values, sizeof(RNScalar), descriptor.nvalues, fp);
    }
  }

  // Write correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    int feature_index0 = correspondence->Feature(0)->index;
    int feature_index1 = correspondence->Feature(1)->index;
    fwrite(&feature_index0, sizeof(int), 1, fp);
    fwrite(&feature_index1, sizeof(int), 1, fp);
    fwrite(&correspondence->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 13; j++) fwrite(&dummy, sizeof(int), 1, fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Sift Feature Input 
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadSiftFile(const char *filename)
{
  // Parse filename
  char run_buffer[4096];
  strncpy(run_buffer, filename, 4096);
  char *file_name_start = run_buffer;
  char *run_name_end = strrchr(run_buffer, '/');
  if (run_name_end) {
    file_name_start = run_name_end + 1;
    *run_name_end = '\0';
  }
  else {
    fprintf(stderr, "Unable to parse run directory from sift filename %s\n", filename); 
    return 0; 
  }

  // Get run name
  printf("%s\n", run_buffer);
  char *run_name = strrchr(run_buffer, '/');
  if (run_name) run_name++;
  else run_name = run_buffer;
  if (!run_name) {
    fprintf(stderr, "Bad run name in sift filename %s\n", filename); 
    return 0; 
  }

  // Parse segment, panorama, and image indicies from filename
  int segment_index, panorama_index, image_index;
  if (sscanf(file_name_start, "%d_%d_%d_Sift.sft",  
    &segment_index, &panorama_index, &image_index) != (unsigned int) 3) {
    fprintf(stderr, "Unable to parse image info from sift filename %s\n", filename);
    return 0;
  }

  // Get run
  GSVRun *run = scene->Run(run_name);
  if (!run) { 
    fprintf(stderr, "Bad run %s in sift filename %s\n", run_name, filename);
    return 0;
  }

  // Get segment index
  if ((segment_index < 0) || (segment_index >= run->NSegments())) { 
    fprintf(stderr, "Bad segment index %d in sift filename %s\n", segment_index, filename);
    return 0;
  }

  // Get segment
  GSVSegment *segment = run->Segment(segment_index);
  if (!segment) {
    fprintf(stderr, "Bad segment %d in sift filename %s\n", segment_index, filename);
    return 0;
  }

  // Get panorama index
  if ((panorama_index < 0) || (panorama_index >= segment->NPanoramas())) { 
    fprintf(stderr, "Bad panorama index %d in sift filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get panorama
  GSVPanorama *panorama = segment->Panorama(panorama_index);
  if (!panorama) {
    fprintf(stderr, "Bad panorama %d in sift filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get image index
  if ((image_index < 0) || (image_index >= panorama->NImages())) {
    fprintf(stderr, "Bad image index %d in sift filename %s\n", image_index, filename);
    return 0;
  }

  // Get image
  GSVImage *image = panorama->Image(image_index);
  if (!image) {
    fprintf(stderr, "Bad image %d in sift filename %s\n", image_index, filename);
    return 0;
  }

  // Read sift file
  return ReadSiftFile(image, filename);
}




int GSVPoseOptimization::
ReadSiftFile(GSVImage *image, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nfeatures, descriptor_size;
  if (fscanf(fp, "%d%d", &nfeatures, &descriptor_size) != (unsigned int) 2) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read features
  double px, py, scale, orientation;
  for (int i = 0; i < nfeatures; i++) {
    // Read feature info
    if (fscanf(fp, "%lf%lf%lf%lf", &py, &px, &scale, &orientation) != (unsigned int) 4) {
      fprintf(stderr, "Unable to read feature %d from %s\n", i, filename);
      return 0;
    }

    // Apply scaling
    px /= sift_image_scale;
    py /= sift_image_scale;
    scale /= sift_image_scale;

    // Invert y coordinate
    py = image->Height() - py;

    // Create feature
    RNScalar score = scale;
    R2Point position(px, py);
    R2Vector direction(1, 0);
    direction.Rotate(orientation);
    GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, image, position, direction, scale, score);
    if (!feature) {
      fprintf(stderr, "Unable to create feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Read feature descriptor
    if (descriptor_size > 0) {
      GSVDescriptor& descriptor = feature->descriptor;
      descriptor.descriptor_type = GSV_SIFT_DESCRIPTOR_TYPE;
      descriptor.nvalues = descriptor_size;
      descriptor.values = new RNScalar [ descriptor_size ];
      for (int j = 0; j < descriptor_size; j++) {
        if (fscanf(fp, "%lf", &descriptor.values[j]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to read descriptor value %d for feature %d in %s\n", j, i, filename);
          return 0;
        }
      }
    }

    // Insert feature
    InsertFeature(feature);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Line Feature Input 
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadLineFile(const char *filename)
{
  // Parse filename
  char run_buffer[4096];
  strncpy(run_buffer, filename, 4096);
  char *file_name_start = run_buffer;
  char *run_name_end = strrchr(run_buffer, '/');
  if (run_name_end) {
    file_name_start = run_name_end + 1;
    *run_name_end = '\0';
  }
  else {
    fprintf(stderr, "Unable to parse run directory from line filename %s\n", filename); 
    return 0; 
  }

  // Get run name
  printf("%s\n", run_buffer);
  char *run_name = strrchr(run_buffer, '/');
  if (run_name) run_name++;
  else run_name = run_buffer;
  if (!run_name) {
    fprintf(stderr, "Bad run name in line filename %s\n", filename); 
    return 0; 
  }

  // Parse segment, panorama, and image indicies from filename
  int segment_index, panorama_index, image_index;
  if (sscanf(file_name_start, "%d_%d_%d_Line.lin",  
    &segment_index, &panorama_index, &image_index) != (unsigned int) 3) {
    fprintf(stderr, "Unable to parse image info from line filename %s\n", filename);
    return 0;
  }

  // Get run
  GSVRun *run = scene->Run(run_name);
  if (!run) { 
    fprintf(stderr, "Bad run %s in line filename %s\n", run_name, filename);
    return 0;
  }

  // Get segment index
  if ((segment_index < 0) || (segment_index >= run->NSegments())) { 
    fprintf(stderr, "Bad segment index %d in line filename %s\n", segment_index, filename);
    return 0;
  }

  // Get segment
  GSVSegment *segment = run->Segment(segment_index);
  if (!segment) {
    fprintf(stderr, "Bad segment %d in line filename %s\n", segment_index, filename);
    return 0;
  }

  // Get panorama index
  if ((panorama_index < 0) || (panorama_index >= segment->NPanoramas())) { 
    fprintf(stderr, "Bad panorama index %d in line filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get panorama
  GSVPanorama *panorama = segment->Panorama(panorama_index);
  if (!panorama) {
    fprintf(stderr, "Bad panorama %d in line filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get image index
  if ((image_index < 0) || (image_index >= panorama->NImages())) {
    fprintf(stderr, "Bad image index %d in line filename %s\n", image_index, filename);
    return 0;
  }

  // Get image
  GSVImage *image = panorama->Image(image_index);
  if (!image) {
    fprintf(stderr, "Bad image %d in line filename %s\n", image_index, filename);
    return 0;
  }

  // Read line file
  return ReadLineFile(image, filename);
}




int GSVPoseOptimization::
ReadLineFile(GSVImage *image, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nfeatures, descriptor_size;
  if (fscanf(fp, "%d%d", &nfeatures, &descriptor_size) != (unsigned int) 2) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read features
  double px1, py1, px2, py2, length, score;
  for (int i = 0; i < nfeatures; i++) {
    // Read feature info
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &px1, &py1, &px2, &py2, &length, &score) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read feature %d from %s\n", i, filename);
      return 0;
    }

    // Read feature descriptor
    GSVDescriptor descriptor;
    if (descriptor_size > 0) {
      descriptor.descriptor_type = GSV_LINE_DESCRIPTOR_TYPE;
      descriptor.nvalues = descriptor_size;
      descriptor.values = new RNScalar [ descriptor_size ];
      for (int j = 0; j < descriptor_size; j++) {
        if (fscanf(fp, "%lf", &descriptor.values[j]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to read descriptor value %d for feature %d in %s\n", j, i, filename);
          return 0;
        }
      }
    }

    // Get useful variables
    R2Point p1(px1, py1);
    R2Point p2(px2, py2);
    R2Vector direction = p2 - p1;
    direction.Normalize();
    RNScalar scale = length;

    // Create feature1
    GSVFeature *feature1 = new GSVFeature(GSV_IMAGE_LINE_FEATURE_TYPE, image, p1, direction, scale, score);
    if (!feature1) {
      fprintf(stderr, "Unable to create line feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Create feature2
    GSVFeature *feature2 = new GSVFeature(GSV_IMAGE_LINE_FEATURE_TYPE, image, p2, -direction, scale, score);
    if (!feature2) {
      fprintf(stderr, "Unable to create line feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Assign feature descriptor
    feature1->descriptor = descriptor;
    feature2->descriptor = descriptor;
      
    // Insert features
    InsertFeature(feature1);
    InsertFeature(feature2);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Transformations Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadTransformationsFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nlasers, ncameras, nvertices, dummy;
  if (fscanf(fp, "%d%d%d%d%d%d%d%d\n", &nlasers, &ncameras, &nvertices, 
    &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 8) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read laser variables 
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    laser_data->translation.Reset(tx, ty, tz);
    // laser_data->rotation.Reset(rx, ry, rz);
  }

  // Read camera variables 
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    camera_data->translation.Reset(tx, ty, tz);
    // camera_data->rotation.Reset(rx, ry, rz);
  }

  // Read path vertex variables
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    vertex->translation.Reset(tx, ty, tz);
    vertex->rotation.Reset(tx, ty, tz);
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteTransformationsFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "%d %d %d 0 0 0 0 0\n", lasers.NEntries(), cameras.NEntries(), vertices.NEntries());

  // Write laser variables to file
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;
    const R3Vector& translation = laser_data->translation;
    const R3Vector& rotation = R3zero_vector; // laser_data->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Write camera variables to file
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;
    const R3Vector& translation = camera_data->translation;
    const R3Vector& rotation = R3zero_vector; // camera_data->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Write path vertex variables to file
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);
    R3Vector translation = vertex->translation;
    R3Vector rotation = vertex->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Update utility functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ApplyOptimizedTransformationsToScene(void)
{
  // Check scene
  if (!scene) return 0;

#if 1
  RNAbort("Not implemented");
#else
  // Set new pose in every scanline
  for (int i = 0; i < scanlines.NEntries(); i++) {
    GSVScanline *scanline = scanlines.Kth(i);
    xxx
    scanline->SetPose(pose);
  }

  // Set new pose in every image
  for (int i = 0; i < images.NEntries(); i++) {
    GSVImage *image = images.Kth(i);
    xxx
    image->SetPose(pose0, pose1);
  }
#endif

  // Return success
  return 1;
}



int GSVPoseOptimization::
UpdatePixelPositionsFromCorrespondenceSolutions(void) 
{
  // Update scan positions for image features
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    if (!feature0 || !feature1) continue;
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      GSVImage *image1 = feature1->image;
      if (!scanline0 || !image1) continue;
      R3Affine transformation0 = OptimizedTransformation(scanline0);
      R3Affine transformation1 = OptimizedTransformation(image1);
      R3Point point0 = feature0->scan_position;
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      R3Ray ray1a = ray1;
      point0.Transform(transformation0);
      ray1.Transform(transformation1);
      RNScalar t = ray1.T(point0);
      feature1->scan_position = ray1a.Point(t);
      feature1->image_t = t;
    }
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      GSVScanline *scanline1 = feature1->scanline;
      if (!image0 || !scanline1) continue;
      R3Affine transformation0 = OptimizedTransformation(image0);
      R3Affine transformation1 = OptimizedTransformation(scanline1);
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Ray ray0a = ray0;
      R3Point point1 = feature1->scan_position;
      ray0.Transform(transformation0);
      point1.Transform(transformation1);
      RNScalar t = ray0.T(point1);
      feature0->scan_position = ray0a.Point(t);
      feature0->image_t = t;
    }
    else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      GSVImage *image1 = feature1->image;
      if (!image0 || !image1) continue;
      R3Affine transformation0 = OptimizedTransformation(image0);
      R3Affine transformation1 = OptimizedTransformation(image1);
      R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
      R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
      R3Ray ray0a = ray0;
      R3Ray ray1a = ray1;
      ray0.Transform(transformation0);
      ray1.Transform(transformation1);
      const R3Vector& v0 = ray0.Vector();
      const R3Vector& v1 = ray1.Vector();
      RNScalar v0v0 = 1.0;  // v0.Dot(v0);
      RNScalar v1v1 = 1.0;  // v1.Dot(v1);
      RNScalar v0v1 = v0.Dot(v1);
      RNScalar denom = v0v1*v0v1 - v0v0*v1v1;
      if (RNIsNotZero(denom)) {
	R3Vector p0 = ray0.Start().Vector();
	R3Vector p1 = ray1.Start().Vector();
	RNScalar p0v0 = v0.Dot(p0);
	RNScalar p1v1 = v1.Dot(p1);
	RNScalar p0v1 = v1.Dot(p0);
	RNScalar p1v0 = v0.Dot(p1);
	RNScalar t0 = (v0v1*p1v1 + v1v1*p0v0 - v0v1*p0v1 - v1v1*p1v0) / denom;
	RNScalar t1 = (v0v1*p0v0 + v0v0*p1v1 - v0v1*p1v0 - v0v0*p0v1) / denom;
	feature0->scan_position = ray0a.Start() + t0 * ray0a.Vector();
	feature1->scan_position = ray1a.Start() + t1 * ray1a.Vector();
        feature0->image_t = t0;
        feature1->image_t = t1;
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
UpdatePixelPositionsFromRayIntersections(void) 
{
  // Parameters
  RNLength max_distance = 100;
  RNScalar max_normal_angle = RN_PI/6.0;
  RNScalar max_curvature = 0;
  RNScalar min_cos_normal_angle = cos(max_normal_angle);

  // Consider features scan-by-scan so that fetch mesh once for each scan
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;
        printf("UPR %d %d %d\n", ir, is, ia);

        // Process all features in segment
        for (int k = 0; k < NFeatures(); k++) {
          GSVFeature *feature = Feature(k);
          GSVImage *image = feature->image;
          if (!image) continue;
          GSVPanorama *panorama = image->Panorama();
          if (!panorama) continue;
          if (panorama->Segment() != segment) continue;

          // Get ray through pixel
          R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
          R3Affine transformation = OptimizedTransformation(image);
          ray.Transform(transformation);

          // Find ray-mesh intersection
          R3MeshIntersection intersection;
          if (!mesh->Intersection(ray, &intersection, 0, max_distance)) continue;

          // Check if closest intersection
          if ((feature->image_t > 0) && (intersection.t > feature->image_t)) continue;
          
          // Get closest mesh vertex
          GSVMeshVertex *vertex = (GSVMeshVertex *) intersection.vertex;
          if (!vertex && intersection.edge) vertex = mesh->VertexOnEdge(intersection.edge);
          if (!vertex && intersection.face) vertex = mesh->VertexOnFace(intersection.face);
          if (!vertex) continue;
          
          // Check normal angle
          if (min_cos_normal_angle > 0) {
            RNScalar dot = -1.0 * ray.Vector().Dot(mesh->VertexNormal(vertex));
            if (dot < min_cos_normal_angle) continue;
          }
              
          // Check curvature
          if (max_curvature > 0) {
            RNScalar curvature = mesh->VertexCurvature(vertex);
            if (curvature > max_curvature) continue;
          }

          // Update feature info
          // GSVScanline *hit_scanline = mesh->VertexScanline(vertex);
          // int hit_point_index = mesh->VertexPointIndex(vertex);
          feature->scan_position = intersection.point;
          feature->scan_normal = mesh->VertexNormal(vertex);
          feature->image_t = intersection.t;

          // Determine scan direction
          if ((feature->image_direction.X() != RN_UNKNOWN) && (feature->image_direction.Y() != RN_UNKNOWN)) {
            // Reset feature info (require ray intersections at both endpoints)
            feature->image_t = 0;
            
            // Get ray intersection at second endpoint
            R2Point image_position2 = feature->image_position + feature->image_scale * feature->image_direction;
            R3Ray ray2 = feature->image->RayThroughUndistortedPosition(image_position2);
            ray2.Transform(transformation);
            R3MeshIntersection intersection2;
            if (!mesh->Intersection(ray2, &intersection2, 0, max_distance)) continue;
            
            // Get closest mesh vertex
            GSVMeshVertex *vertex2 = (GSVMeshVertex *) intersection2.vertex;
            if (!vertex2 && intersection2.edge) vertex2 = mesh->VertexOnEdge(intersection2.edge);
            if (!vertex2 && intersection2.face) vertex2 = mesh->VertexOnFace(intersection2.face);
            if (!vertex2) continue;
            
            // Check normal angle
            if (min_cos_normal_angle > 0) {
              RNScalar dot2 = -1.0 * ray2.Vector().Dot(mesh->VertexNormal(vertex2));
              if (dot2 < min_cos_normal_angle) continue;
            }
            
            // Check curvature
            if (max_curvature > 0) {
              RNScalar curvature2 = mesh->VertexCurvature(vertex2);
              if (curvature2 > max_curvature) continue;
            }
            
            // Update feature info
            feature->scan_position = intersection.point;
            feature->scan_normal = mesh->VertexNormal(vertex);
            feature->scan_direction = intersection2.point - intersection.point;
            feature->scan_scale = feature->scan_direction.Length();
            if (RNIsNotZero(feature->scan_scale)) feature->scan_direction /= feature->scan_scale;
            feature->image_t = intersection.t;
          }
        }

        // Delete mesh
        delete mesh;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////
// Polynomial equation stuff
////////////////////////////

void GSVPoseOptimization::
AddInertiaEquations(RNPolynomialSystemOfEquations *system)
{
  // Insert equations for laser variables
  if (laser_nv > 0) {
    for (int i = 0; i < lasers.NEntries(); i++) {
      for (int j = 0; j < LASER_NV; j++) {
        if (laser_v[j] == -1) continue;
        RNPolynomial *equation = new RNPolynomial();
        AddLaserVariableToEquation(equation, i, j, laser_inertia_weights[j]);
        system->InsertPolynomial(equation);
      }
    }
  }

  // Insert equations for camera variables
  if (camera_nv > 0) {
    for (int i = 0; i < cameras.NEntries(); i++) {
      for (int j = 0; j < CAMERA_NV; j++) {
        if (camera_v[j] == -1) continue;
        RNPolynomial *equation = new RNPolynomial();
        AddCameraVariableToEquation(equation, i, j, camera_inertia_weights[j]);
        system->InsertPolynomial(equation);
      }
    }
  }

  // Insert equations for path variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      for (int j = 0; j < VERTEX_NV; j++) {
        if (vertex_v[j] == -1) continue;
        RNPolynomial *equation = new RNPolynomial();
        AddPathVertexVariableToEquation(equation, i, j, vertex_inertia_weights[j]);
        system->InsertPolynomial(equation);
      }
    }
  }

  // Insert equations for pixel variables
  if (pixel_nv > 0) {
    for (int i = 0; i < pixels.NEntries(); i++) {
      for (int j = 0; j < PIXEL_NV; j++) {
        if (pixel_v[j] == -1) continue;
        RNPolynomial *equation = new RNPolynomial();
        AddPixelVariableToEquation(equation, i, j, pixel_inertia_weights[j]);
        system->InsertPolynomial(equation);
      }
    }
  }
}



void GSVPoseOptimization::
AddRigidityEquations(RNPolynomialSystemOfEquations *system, RNScalar rigidity)
{
  // Check if optimizing path variables
  if (vertex_nv == 0) return;

  // Insert equations to maintain the shape of the path
  RNScalar rigidity_scale_factor = 100 * pow(0.5 + rigidity, 5);
  for (int is = 0; is < segments.NEntries(); is++) {
    GSVSegment *segment = segments.Kth(is);
    SegmentData *segment_data = (SegmentData *) segment->Data();
    if (!segment_data) continue;
    GSVPath *path = segment_data->path;
    if (!path) continue;
    for (int i0 = 0; i0 < path->NVertices(); i0++) {
      GSVPathVertex *vertex0 = path->Vertex(i0);
      const R3Point& position0 = path->VertexPosition(i0);
      for (int k = 1; k <= intra_segment_rigidity_radius; k++) {
        int i1 = i0 + k;
        if (i1 >= path->NVertices()) continue;
        GSVPathVertex *vertex1 = path->Vertex(i1);
        const R3Point& position1 = path->VertexPosition(i1);
        RNScalar s = exp(-2 * k * k / (intra_segment_rigidity_sigma_squared));
        RNScalar w = s * rigidity_scale_factor * intra_segment_rigidity_weight;
        AddPathRigidityEquations(system, vertex0->index, vertex1->index, position0, position1, w);
        AddPathRigidityEquations(system, vertex1->index, vertex0->index, position1, position0, w);
      }
    }
  }
}



void GSVPoseOptimization::
AddPathRigidityEquations(RNPolynomialSystemOfEquations *system,
  int vertex_index0, int vertex_index1, 
  const R3Point& v0, const R3Point& v1, RNScalar w)
{
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |

  // (R0(d) + v0 + t0) - (v1 + t1) = 0
  // (R0(d) - d + t0 - t1 = 0

  // X:      d.X() + -rz0*d.Y() +  ry0*d.Z() + -d.X() + tx0 - s1tx = 0
  // Y:  rz0*d.X() +      d.Y() + -rx0*d.Z() + -d.Y() + ty0 - s1ty = 0
  // Z: -ry0*d.X() +  rx0*d.Y() +      d.Z() + -d.Z() + tz0 - s1tz = 0

  // X:      -rz0*d.Y() +   ry0*d.Z() + tx0 - s1tx = 0
  // Y:       rz0*d.X() +  -rx0*d.Z() + ty0 - s1ty = 0
  // Z:      -ry0*d.X() +   rx0*d.Y() + tz0 - s1tz = 0

  // Check if optimizing path variables
  if (vertex_nv == 0) return;

  // Get convenient variables
  RNPolynomial *equation;
  R3Vector d = v1 - v0;

  // Align X coordinates
  equation = new RNPolynomial();
  AddPathVertexVariableToEquation(equation, vertex_index0, TX,  1.0); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RY,  d.Z()); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RZ, -d.Y()); 
  AddPathVertexVariableToEquation(equation, vertex_index1, TX, -1.0); 
  equation->Multiply(w);                              
  system->InsertPolynomial(equation);                 
                                                      
  // Align Y coordinates                              
  equation = new RNPolynomial();                      
  AddPathVertexVariableToEquation(equation, vertex_index0, TY,  1.0); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RX, -d.Z()); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RZ,  d.X()); 
  AddPathVertexVariableToEquation(equation, vertex_index1, TY, -1.0); 
  equation->Multiply(w);                              
  system->InsertPolynomial(equation);                 
                                                      
  // Align Z coordinates                              
  equation = new RNPolynomial();                      
  AddPathVertexVariableToEquation(equation, vertex_index0, TZ,  1.0);      
  AddPathVertexVariableToEquation(equation, vertex_index0, RX,  d.Y());
  AddPathVertexVariableToEquation(equation, vertex_index0, RY, -d.X());
  AddPathVertexVariableToEquation(equation, vertex_index1, TZ, -1.0);      
  equation->Multiply(w);
  system->InsertPolynomial(equation);
}



void GSVPoseOptimization::
AddCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
  RNBoolean include_scan_scan_correspondences, 
  RNBoolean include_image_image_correspondences, 
  RNBoolean include_scan_image_correspondences)
{
  // Encourage correspondences 
  for (int c = 0; c < NCorrespondences(); c++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(c);
    AddCorrespondenceEquations(system, correspondence, 
      include_scan_scan_correspondences, 
      include_image_image_correspondences, 
      include_scan_image_correspondences);
  }
}



void GSVPoseOptimization::
AddCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVFeatureCorrespondence *correspondence,
  RNBoolean include_scan_scan_correspondences, 
  RNBoolean include_image_image_correspondences, 
  RNBoolean include_scan_image_correspondences)
{
  // Get features
  GSVFeature *feature0 = correspondence->Feature(0);
  GSVFeature *feature1 = correspondence->Feature(1);
  if (!feature0 || !feature1) return;

  // SCAN <-> SCAN
  if (include_scan_scan_correspondences) {
    // SCANPOINT <-> SCANPOINT
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      RNScalar w = scan_point_scan_point_correspondence_weight * correspondence->score;
      AddPointPointCorrespondenceEquations(system, scanline0, point0, scanline1, point1, w);
    }

    // SCANPOINT <-> SCANPLANE
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      R3Vector normal1 = feature1->scan_normal;
      RNScalar w = scan_point_scan_plane_correspondence_weight * correspondence->score;
      AddPointPlaneCorrespondenceEquations(system, scanline0, point0, scanline1, point1, normal1, w);
    }

    // SCANPLANE <-> SCANPOINT
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      R3Vector normal0 = feature0->scan_normal;
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      RNScalar w = scan_point_scan_plane_correspondence_weight * correspondence->score;
      AddPointPlaneCorrespondenceEquations(system, scanline1, point1, scanline0, point0, normal0, w);
    }

    // SCANPLANE <-> SCANPLANE
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      R3Vector normal0 = feature0->scan_normal;
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      R3Vector normal1 = feature1->scan_normal;
      RNScalar w = scan_point_scan_plane_correspondence_weight * correspondence->score;
      AddPointPlaneCorrespondenceEquations(system, scanline1, point1, scanline0, point0, normal0, w);
      AddPointPlaneCorrespondenceEquations(system, scanline0, point0, scanline1, point1, normal1, w);
    }
  }

  // IMAGE <-> IMAGE
  if (include_image_image_correspondences) {
    // IMAGEPOINT -> IMAGEPOINT
    if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      R2Point point0 = feature0->image_position;
      int pixel_index0 = feature0->pixel_index;
      R3Ray ray0 = image0->RayThroughUndistortedPosition(point0);
      GSVImage *image1 = feature1->image;
      R2Point point1 = feature1->image_position;
      int pixel_index1 = feature1->pixel_index;
      R3Ray ray1 = image1->RayThroughUndistortedPosition(point1);
      RNScalar w = image_point_image_point_correspondence_weight * correspondence->score;
      AddPixelPixelCorrespondenceEquations(system, image0, point0, pixel_index0, image1, point1, pixel_index1, w);
      AddPixelRayCorrespondenceEquations(system, image0, point0, pixel_index0, image1, ray1, w);
      AddPixelRayCorrespondenceEquations(system, image1, point1, pixel_index1, image0, ray0, w);
    }
  }

  // SCAN <-> IMAGE
  if (include_scan_image_correspondences) {
    // SCANPOINT <-> IMAGEPOINT
    if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      GSVImage *image1 = feature1->image;
      R2Point point1 = feature1->image_position;
      int pixel_index1 = feature1->pixel_index;
      R3Ray ray1 = image1->RayThroughUndistortedPosition(point1);
      RNScalar w = scan_point_image_point_correspondence_weight * correspondence->score;
      AddPointRayCorrespondenceEquations(system, scanline0, point0, image1, ray1, w);
      AddPixelPointCorrespondenceEquations(system, image1, point1, pixel_index1, scanline0, point0, w);
    }

    // SCANPLANE <-> IMAGEPOINT
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVScanline *scanline0 = feature0->scanline;
      R3Point point0 = feature0->scan_position;
      R3Vector normal0 = feature0->scan_normal;
      GSVImage *image1 = feature1->image;
      R2Point point1 = feature1->image_position;
      int pixel_index1 = feature1->pixel_index;
      RNScalar w = image_point_scan_plane_correspondence_weight * correspondence->score;
      AddPixelPlaneCorrespondenceEquations(system, image1, point1, pixel_index1, scanline0, point0, normal0, w);
    }

    // IMAGEPOINT <-> SCANPOINT
    if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      R2Point point0 = feature0->image_position;
      int pixel_index0 = feature0->pixel_index;
      R3Ray ray0 = image0->RayThroughUndistortedPosition(point0);
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      RNScalar w = scan_point_image_point_correspondence_weight * correspondence->score;
      AddPointRayCorrespondenceEquations(system, scanline1, point1, image0, ray0, w);
      AddPixelPointCorrespondenceEquations(system, image0, point0, pixel_index0, scanline1, point1, w);
    }

    // IMAGEPOINT <-> SCANPLANE 
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      GSVImage *image0 = feature0->image;
      R2Point point0 = feature0->image_position;
      int pixel_index0 = feature0->pixel_index;
      R3Ray ray0 = image0->RayThroughUndistortedPosition(point0);
      GSVScanline *scanline1 = feature1->scanline;
      R3Point point1 = feature1->scan_position;
      R3Vector normal1 = feature1->scan_normal;
      RNScalar w = image_point_scan_plane_correspondence_weight * correspondence->score;
      AddPixelPlaneCorrespondenceEquations(system, image0, point0, pixel_index0, scanline1, point1, normal1, w);
    }
  }
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVScanline *scanline, const R3Point& point,
  RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz)
{
  // Get useful variables
  if (!scanline) return 0;
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return 0;
  int scanline_index = scanline_data->index;
  if (scanline_index < 0) return 0;
  GSVScan *scan = scanline->Scan();
  if (!scan) return 0;
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get path pose in world coordinates
  const GSVPose& path_pose = path->Pose(scanline_data->path_parameter);
  const R3Point& path_viewpoint = path_pose.Viewpoint();

  // Get vector from path viewpoint to point
  RNPolynomial dx, dy, dz;
  dx += point.X() - path_viewpoint.X();
  dy += point.Y() - path_viewpoint.Y();
  dz += point.Z() - path_viewpoint.Z();

  // Add laser transformation
  if (laser_nv > 0) {
    // Get more convenient variables
    GSVLaser *laser = scanline->Scan()->Laser();
    if (!laser) return 0;
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) return 0;
    int laser_index = laser_data->index;
    if (laser_index < 0) return 0;
    const R3Vector path_right = path_pose.Right();
    const R3Vector path_towards = path_pose.Towards();
    const R3Vector path_up = path_pose.Up();

    // Add laser transformation (in path coordinate system)
    RNPolynomial ltx, lty, ltz;
    AddLaserVariableToEquation(&ltx, laser_index, TX); 
    AddLaserVariableToEquation(&lty, laser_index, TY); 
    AddLaserVariableToEquation(&ltz, laser_index, TZ); 
    dx += ltx*path_right.X(); dx += lty*path_towards.X(); dx += ltz*path_up.X();
    dy += ltx*path_right.Y(); dy += lty*path_towards.Y(); dy += ltz*path_up.Y();
    dz += ltx*path_right.Z(); dz += lty*path_towards.Z(); dz += ltz*path_up.Z();
  }

  // Get path transformation variables 
  RNPolynomial stx, sty, stz, srx, sry, srz;
  AddScanlineVariableToEquation(&stx, scanline_index, TX); 
  AddScanlineVariableToEquation(&sty, scanline_index, TY); 
  AddScanlineVariableToEquation(&stz, scanline_index, TZ); 
  AddScanlineVariableToEquation(&srx, scanline_index, RX); 
  AddScanlineVariableToEquation(&sry, scanline_index, RY); 
  AddScanlineVariableToEquation(&srz, scanline_index, RZ); 

  // Get coordinates for transformed point
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |
  px = stx;
  py = sty;
  pz = stz;
  px += path_viewpoint.X();
  py += path_viewpoint.Y();
  pz += path_viewpoint.Z();
  px +=     dx;  px -= srz*dy;  px += sry*dz;
  py += srz*dx;  py +=     dy;  py -= srx*dz;
  pz -= sry*dx;  pz += srx*dy;  pz +=     dz;

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVImage *image, const R3Point& point,
  RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz)
{
  // Get useful variables
  if (!image) return 0;
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return 0;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get vector from path viewpoint to point
  RNPolynomial dx, dy, dz;
  const GSVPose& path_pose = path->Pose(image_data->path_parameter);
  const R3Point& path_viewpoint = path_pose.Viewpoint();
  dx += point.X() - path_viewpoint.X();
  dy += point.Y() - path_viewpoint.Y();
  dz += point.Z() - path_viewpoint.Z();

  // Add camera transformation
  if (camera_nv > 0) {
    // Get more convenient variables
    GSVCamera *camera = image->Tapestry()->Camera();
    if (!camera) return 0;
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) return 0;
    int camera_index = camera_data->index;
    if (camera_index < 0) return 0;
    const R3Vector path_right = path_pose.Right();
    const R3Vector path_towards = path_pose.Towards();
    const R3Vector path_up = path_pose.Up();

    // Add camera transformation (in path coordinate system)
    RNPolynomial ctx, cty, ctz;
    AddCameraVariableToEquation(&ctx, camera_index, TX); 
    AddCameraVariableToEquation(&cty, camera_index, TY); 
    AddCameraVariableToEquation(&ctz, camera_index, TZ); 
    dx += ctx*path_right.X(); dx += cty*path_towards.X(); dx += ctz*path_up.X();
    dy += ctx*path_right.Y(); dy += cty*path_towards.Y(); dy += ctz*path_up.Y();
    dz += ctx*path_right.Z(); dz += cty*path_towards.Z(); dz += ctz*path_up.Z();
  }

  // Get path transformation variables 
  RNPolynomial stx, sty, stz, srx, sry, srz;
  AddImageVariableToEquation(&stx, image_index, TX); 
  AddImageVariableToEquation(&sty, image_index, TY); 
  AddImageVariableToEquation(&stz, image_index, TZ); 
  AddImageVariableToEquation(&srx, image_index, RX); 
  AddImageVariableToEquation(&sry, image_index, RY); 
  AddImageVariableToEquation(&srz, image_index, RZ); 

  // Get coordinates for transformed point
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |
  px = stx; 
  py = sty; 
  pz = stz; 
  px += path_viewpoint.X();
  py += path_viewpoint.Y();
  pz += path_viewpoint.Z();
  px +=     dx;  px -= srz*dy;  px += sry*dz;
  py += srz*dx;  py +=     dy;  py -= srx*dz;
  pz -= sry*dx;  pz += srx*dy;  pz +=     dz;

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVImage *image, const R2Point& point, int pixel_index,
  RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz)
{
  // Get useful variables
  if (!image) return 0;
  if (pixel_index < 0) return 0;
  R3Ray ray = image->RayThroughUndistortedPosition(point);
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return 0;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get ray start offset from path point
  RNPolynomial ox, oy, oz;
  const R3Point& v = ray.Start();
  const GSVPose& path_pose = path->Pose(image_data->path_parameter);
  const R3Point& path_viewpoint = path_pose.Viewpoint();
  ox += v.X() - path_viewpoint.X();
  oy += v.Y() - path_viewpoint.Y();
  oz += v.Z() - path_viewpoint.Z();

  // Add camera transformation
  if (camera_nv > 0) {
    // Get more convenient variables
    GSVCamera *camera = image->Tapestry()->Camera();
    if (!camera) return 0;
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) return 0;
    int camera_index = camera_data->index;
    if (camera_index < 0) return 0;
    const R3Vector path_right = path_pose.Right();
    const R3Vector path_towards = path_pose.Towards();
    const R3Vector path_up = path_pose.Up();
    
    // Add camera transformation (in path coordinate system)
    RNPolynomial ctx, cty, ctz;
    AddCameraVariableToEquation(&ctx, camera_index, TX); 
    AddCameraVariableToEquation(&cty, camera_index, TY); 
    AddCameraVariableToEquation(&ctz, camera_index, TZ); 
    ox = ctx*path_right.X(); ox += cty*path_towards.X(); ox += ctz*path_up.X();
    oy = ctx*path_right.Y(); oy += cty*path_towards.Y(); oy += ctz*path_up.Y();
    oz = ctx*path_right.Z(); oz += cty*path_towards.Z(); oz += ctz*path_up.Z();
  }

  // Get path transformation variables 
  RNPolynomial stx, sty, stz, srx, sry, srz;
  AddImageVariableToEquation(&stx, image_index, TX); 
  AddImageVariableToEquation(&sty, image_index, TY); 
  AddImageVariableToEquation(&stz, image_index, TZ); 
  AddImageVariableToEquation(&srx, image_index, RX); 
  AddImageVariableToEquation(&sry, image_index, RY); 
  AddImageVariableToEquation(&srz, image_index, RZ); 

  // Get pixel depth 
  RNPolynomial t;
  AddPixelVariableToEquation(&t, pixel_index, T, 1.0); 

  // Get transformed ray start
  RNPolynomial vx = stx;
  RNPolynomial vy = sty;
  RNPolynomial vz = stz;
  vx += path_viewpoint.X();
  vy += path_viewpoint.Y();
  vz += path_viewpoint.Z();
  vx +=     ox;  vx -= srz*oy;  vx += sry*oz;   
  vy += srz*ox;  vy +=     oy;  vy -= srx*oz;
  vz -= sry*ox;  vz += srx*oy;  vz +=   oz;

  // Get transformed ray vector
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |
  const R3Vector& d = ray.Vector();
  RNPolynomial dx =      d.X() - srz*d.Y() + sry*d.Z();
  RNPolynomial dy =  srz*d.X() +     d.Y() - srx*d.Z();
  RNPolynomial dz = -sry*d.X() + srx*d.Y() +     d.Z();

  // Get point position
  px = vx;  px += t*dx;
  py = vy;  py += t*dy;
  pz = vz;  pz += t*dz;

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedVectorCoordinates(GSVScanline *scanline, const R3Vector& d,
  RNPolynomial& dx, RNPolynomial& dy, RNPolynomial& dz)
{
  // sr = rotation of spline path around s
  // vr = rotation around v
  // d' = sr(vr(d))

  // Get useful variables
  if (!scanline) return 0;
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return 0;
  int scanline_index = scanline_data->index;
  if (scanline_index < 0) return 0;

  // Get path rotation
  RNPolynomial srx, sry, srz;
  AddScanlineVariableToEquation(&srx, scanline_index, RX); 
  AddScanlineVariableToEquation(&sry, scanline_index, RY); 
  AddScanlineVariableToEquation(&srz, scanline_index, RZ); 

  // Get rotated vector
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |
  dx =      d.X() - srz*d.Y() + sry*d.Z();
  dy =  srz*d.X() +     d.Y() - srx*d.Z();
  dz = -sry*d.X() + srx*d.Y() +     d.Z();

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedVectorCoordinates(GSVImage *image, const R3Vector& d,
  RNPolynomial& dx, RNPolynomial& dy, RNPolynomial& dz)
{
  // Get useful variables
  if (!image) return 0;
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;

  // Get image transformation variables
  RNPolynomial srx, sry, srz;
  AddImageVariableToEquation(&srx, image_index, RX); 
  AddImageVariableToEquation(&sry, image_index, RY); 
  AddImageVariableToEquation(&srz, image_index, RZ); 

  // Get vector
  dx =      d.X() - srz*d.Y() + sry*d.Z();
  dy =  srz*d.X() +     d.Y() - srx*d.Z();
  dz = -sry*d.X() + srx*d.Y() +     d.Z();

  // Return success
  return 1;
}



void GSVPoseOptimization::
AddPointPointCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVScanline *scanline0, const R3Point& point0, 
  GSVScanline *scanline1, const R3Point& point1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1;
  if (!ComputeTransformedPointCoordinates(scanline0, point0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(scanline1, point1, px1, py1, pz1)) return;
    
  // Add equation representing differences between transformed point coordinates
  RNPolynomial *ex = new RNPolynomial();
  RNPolynomial *ey = new RNPolynomial();
  RNPolynomial *ez = new RNPolynomial();
  *ex += w *(px1 - px0);
  *ey += w *(py1 - py0);
  *ez += w *(pz1 - pz0);
  system->InsertPolynomial(ex);
  system->InsertPolynomial(ey);
  system->InsertPolynomial(ez);
}



void GSVPoseOptimization::
AddPointRayCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
  GSVScanline *scanline0, const R3Point& point0, 
  GSVImage *image1, const R3Ray& ray1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

#if 0
  // Move point1 to most stable position
  R3Point point1 = point0;
  point1.Project(ray1.Line());
  R3Point direction1 = ray1.Vector();

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, dx1, dy1, dz1;
  if (!ComputeTransformedPointCoordinates(scanline0, point0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(image1, point1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(image1, direction1, dx1, dy1, dz1)) return;
    
  // Compute vector from point1 to point0
  RNPolynomial vx = px1 - px0;
  RNPolynomial vy = py1 - py0;
  RNPolynomial vz = pz1 - pz0;
  
  // Add equations for distance between point0 and line1
  // || d % v || = 0;
  RNPolynomial *ex = new RNPolynomial();
  RNPolynomial *ey = new RNPolynomial();
  RNPolynomial *ez = new RNPolynomial();
  *ex = w * (dy1*vz - dz1*vy);
  *ey = w * (dz1*vx - dx1*vz);
  *ez = w * (dx1*vy - dy1*vx);
  system->InsertPolynomial(ex);
  system->InsertPolynomial(ey);
  system->InsertPolynomial(ez);
#else
  // Add equations for orthogonal planes intersecting at ray1
  R3Vector normal1a = R3posz_vector % ray1.Vector(); normal1a.Normalize();
  R3Vector normal1b = normal1a % ray1.Vector(); normal1b.Normalize();
  AddPointPlaneCorrespondenceEquations(system, scanline0, point0, image1, ray1.Start(), normal1a, w);
  AddPointPlaneCorrespondenceEquations(system, scanline0, point0, image1, ray1.Start(), normal1b, w);
#endif
}



void GSVPoseOptimization::
AddPointPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
  GSVScanline *scanline0, const R3Point& point0, 
  GSVScanline *scanline1, const R3Point& point1, const R3Vector& normal1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Find point ont plane1 closest to point0 (adds stability)
  R3Plane plane1(point1, normal1);
  R3Point center1 = point0;
  center1.Project(plane1);

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1;
  if (!ComputeTransformedPointCoordinates(scanline0, point0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(scanline1, center1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(scanline1, normal1, nx1, ny1, nz1)) return;

  // Add equation representing distance from transformed point0 to transformed plane1
  RNPolynomial *equation = new RNPolynomial();
  *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 + (pz0 - pz1)*nz1 );
  system->InsertPolynomial(equation);
}



void GSVPoseOptimization::
AddPointPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
  GSVScanline *scanline0, const R3Point& point0, 
  GSVImage *image1, const R3Point& point1, const R3Vector& normal1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Find point ont plane1 closest to point0 (adds stability)
  R3Plane plane1(point1, normal1);
  R3Point center1 = point0;
  center1.Project(plane1);

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1;
  if (!ComputeTransformedPointCoordinates(scanline0, point0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(image1, center1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(image1, normal1, nx1, ny1, nz1)) return;

  // Compute equation representing distance from transformed point0 to transformed plane1
  RNPolynomial *equation = new RNPolynomial();
  *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 + (pz0 - pz1)*nz1 );
  system->InsertPolynomial(equation);
}



void GSVPoseOptimization::
AddPixelPointCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVImage *image0, const R2Point& point0, int pixel_index0,
  GSVScanline *scanline1, const R3Point& point1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1;
  if (!ComputeTransformedPointCoordinates(image0, point0, pixel_index0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(scanline1, point1, px1, py1, pz1)) return;
    
  // Compute distance between points
  RNPolynomial *ex = new RNPolynomial();
  RNPolynomial *ey = new RNPolynomial();
  RNPolynomial *ez = new RNPolynomial();
  *ex = w * ( px0 - px1 );
  *ey = w * ( py0 - py1 );
  *ez = w * ( pz0 - pz1 );

  // Divide by distance to get approximation to angle
  // RNPolynomial inverse_t0;
  // AddPixelVariableToEquation(&inverse_t0, pixel_index0, T, 1.0, -1.0); 
  // *ex *= inverse_t0;
  // *ey *= inverse_t0;
  // *ez *= inverse_t0;

  // Insert equations
  system->InsertPolynomial(ex);
  system->InsertPolynomial(ey);
  system->InsertPolynomial(ez);
}



void GSVPoseOptimization::
AddPixelRayCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVImage *image0, const R2Point& point0, int pixel_index0, 
  GSVImage *image1, const R3Ray& ray1, 
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

#if 0
  // Move point1 to most stable position
  R3Point point1 = ray1.Start();
  R3Point direction1 = ray1.Vector();
  R3Ray ray0 = image->RayThroughUndistortedPosition(point0);
  RNScalar t0 = PixelVariableValue(pixel_index, T);
  if (t0 > 0) { point1 = ray0.Point(t0); point1.Project(ray1.Line()); }

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, dx1, dy1, dz1;
  if (!ComputeTransformedPointCoordinates(image0, point0, pixel_index0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(image1, point1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(image1, direction1, dx1, dy1, dz1)) return;
    
  // Compute vector from point1 to point0
  RNPolynomial vx = px1 - px0;
  RNPolynomial vy = py1 - py0;
  RNPolynomial vz = pz1 - pz0;
  
  // Add equations for distance between point0 and line1
  // || d % v || = 0;
  RNPolynomial *ex = new RNPolynomial();
  RNPolynomial *ey = new RNPolynomial();
  RNPolynomial *ez = new RNPolynomial();
  *ex = w * (dy1*vz - dz1*vy);
  *ey = w * (dz1*vx - dx1*vz);
  *ez = w * (dx1*vy - dy1*vx);

  // Divide by distance to get approximation to angle
  // RNPolynomial inverse_t0;
  // AddPixelVariableToEquation(&inverse_t0, pixel_index0, T, 1.0, -1.0); 
  // *ex *= inverse_t0;
  // *ey *= inverse_t0;
  // *ez *= inverse_t0;

  // Insert equations
  system->InsertPolynomial(ex);
  system->InsertPolynomial(ey);
  system->InsertPolynomial(ez);
#else
  // Add equations for orthogonal planes intersecting at ray1
  R3Vector normal1a = R3posz_vector % ray1.Vector(); normal1a.Normalize();
  R3Vector normal1b = normal1a % ray1.Vector(); normal1b.Normalize();
  AddPixelPlaneCorrespondenceEquations(system, image0, point0, pixel_index0, image1, ray1.Start(), normal1a, w);
  AddPixelPlaneCorrespondenceEquations(system, image0, point0, pixel_index0, image1, ray1.Start(), normal1b, w);
#endif
}



void GSVPoseOptimization::
AddPixelPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVImage *image0, const R2Point& point0, int pixel_index0,
  GSVScanline *scanline1, const R3Point& point1, const R3Vector& normal1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Find point on plane1 closest to expected position of point0 (adds stability)
  R3Point center1 = point1;
  R3Plane plane1(point1, normal1);
  R3Ray ray0 = image0->RayThroughUndistortedPosition(point0);
  RNScalar t0 = PixelVariableValue(pixel_index0, T);
  if (t0 > 0) { center1 = ray0.Point(t0); center1.Project(plane1); }

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1;
  if (!ComputeTransformedPointCoordinates(image0, point0, pixel_index0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(scanline1, center1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(scanline1, normal1, nx1, ny1, nz1)) return;

  // Add equation representing distance from transformed point1 to transformed plane0
  RNPolynomial *equation = new RNPolynomial();
  *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 + (pz0 - pz1)*nz1 );

  // Divide by distance to get approximation to angle
  // RNPolynomial inverse_t0;
  // AddPixelVariableToEquation(&inverse_t0, pixel_index0, T, 1.0, -1.0); 
  // *equation *= inverse_t0;

  // Insert equation
  system->InsertPolynomial(equation);
}



void GSVPoseOptimization::
AddPixelPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
  GSVImage *image0, const R2Point& point0, int pixel_index0,
  GSVImage *image1, const R3Point& point1, const R3Vector& normal1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Find point on plane1 closest to expected position of point0 (adds stability)
  R3Point center1 = point1;
  R3Plane plane1(point1, normal1);
  R3Ray ray0 = image0->RayThroughUndistortedPosition(point0);
  RNScalar t0 = PixelVariableValue(pixel_index0, T);
  if (t0 > 0) { center1 = ray0.Point(t0); center1.Project(plane1); }

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1;
  if (!ComputeTransformedPointCoordinates(image0, point0, pixel_index0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(image1, center1, px1, py1, pz1)) return;
  if (!ComputeTransformedVectorCoordinates(image1, normal1, nx1, ny1, nz1)) return;

  // Add equation representing distance from transformed point1 to transformed plane0
  RNPolynomial *equation = new RNPolynomial();
  *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 + (pz0 - pz1)*nz1 );

  // Divide by distance to get approximation to angle
  // RNPolynomial inverse_t0;
  // AddPixelVariableToEquation(&inverse_t0, pixel_index0, T, 1.0, -1.0); 
  // *equation *= inverse_t0;

  // Insert equation
  system->InsertPolynomial(equation);
}



void GSVPoseOptimization::
AddPixelPixelCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
  GSVImage *image0, const R2Point& point0, int pixel_index0, 
  GSVImage *image1, const R2Point& point1, int pixel_index1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Get transformed coordinates
  RNPolynomial px0, py0, pz0, px1, py1, pz1;
  if (!ComputeTransformedPointCoordinates(image0, point0, pixel_index0, px0, py0, pz0)) return;
  if (!ComputeTransformedPointCoordinates(image1, point1, pixel_index1, px1, py1, pz1)) return;

  // Add equations for distance between transformed pixels
  RNPolynomial *ex0 = new RNPolynomial();
  RNPolynomial *ey0 = new RNPolynomial();
  RNPolynomial *ez0 = new RNPolynomial();
  RNPolynomial *ex1 = new RNPolynomial();
  RNPolynomial *ey1 = new RNPolynomial();
  RNPolynomial *ez1 = new RNPolynomial();
  RNPolynomial dx = px0 - px1;
  RNPolynomial dy = py0 - py1;
  RNPolynomial dz = pz0 - pz1;
  *ex0 = w * dx;
  *ey0 = w * dy;
  *ez0 = w * dz;
  *ex1 = w * dx;
  *ey1 = w * dy;
  *ez1 = w * dz;

  // Divide by distance to get approximation to angle
  // RNPolynomial inverse_t0, inverse_t1;
  // AddPixelVariableToEquation(&inverse_t0, pixel_index0, T, 1.0, -1.0); 
  // AddPixelVariableToEquation(&inverse_t1, pixel_index1, T, 1.0, -1.0); 
  // *ex0 *= inverse_t0;
  // *ey0 *= inverse_t0;
  // *ez0 *= inverse_t0;
  // *ex1 *= inverse_t1;
  // *ey1 *= inverse_t1;
  // *ez1 *= inverse_t1;

  // Insert equations
  system->InsertPolynomial(ex0);
  system->InsertPolynomial(ey0);
  system->InsertPolynomial(ez0);
  system->InsertPolynomial(ex1);
  system->InsertPolynomial(ey1);
  system->InsertPolynomial(ez1);
}



////////////////////////////////////////////////////////////////////////
// Low level equation functions
////////////////////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddLaserVariableToEquation(RNPolynomial *equation, 
  int laser_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (laser_nv == 0) return;
  assert(equation);
  assert(laser_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < LASER_NV);

  // Add variable term to equation
  int index = LaserVariableIndex(laser_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * LaserVariableValue(laser_index, variable_index));
}



void GSVPoseOptimization::
AddCameraVariableToEquation(RNPolynomial *equation, 
  int camera_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (camera_nv == 0) return;
  assert(equation);
  assert(camera_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < CAMERA_NV);

  // Add variable term to equation
  int index = CameraVariableIndex(camera_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * CameraVariableValue(camera_index, variable_index));
}



void GSVPoseOptimization::
AddPathVertexVariableToEquation(RNPolynomial *equation, 
  int vertex_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (vertex_nv == 0) return;
  assert(equation);
  assert(vertex_index >= 0);
  assert(vertex_index < vertices.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);

  // Check if variable is constrained
  int index = VertexVariableIndex(vertex_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * VertexVariableValue(vertex_index, variable_index));
}



void GSVPoseOptimization::
AddPathVertexVariableToEquation(RNPolynomial *equation, 
  GSVPath *path, RNScalar path_parameter, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (vertex_nv == 0) return;
  assert(equation);
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddScanlineVariableToEquation(RNPolynomial *equation, 
  int scanline_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (vertex_nv == 0) return;
  assert(equation);
  assert(scanline_index >= 0);
  assert(scanline_index < scanlines.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  GSVScanline *scanline = scanlines.Kth(scanline_index);
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return;
  GSVScan *scan = scanline->Scan();
  if (!scan) return;
  GSVSegment *segment = scan->Segment();
  if (!segment) return;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return;
  GSVPath *path = segment_data->path;
  if (!path) return;
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = scanline_data->path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddImageVariableToEquation(RNPolynomial *equation, 
  int image_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (vertex_nv == 0) return;
  assert(equation);
  assert(image_index >= 0);
  assert(image_index < images.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  GSVImage *image = images.Kth(image_index);
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return;
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return;
  GSVSegment *segment = panorama->Segment();
  if (!segment) return;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return;
  GSVPath *path = segment_data->path;
  if (!path) return;
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = image_data->path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddPixelVariableToEquation(RNPolynomial *equation, 
  int pixel_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  if (pixel_nv == 0) return;
  assert(equation);
  assert(pixel_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < PIXEL_NV);

  // Add variable term to equation
  int index = PixelVariableIndex(pixel_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * PixelVariableValue(pixel_index, variable_index));
}



////////////////////////////
// Optimization stuff
////////////////////////////

int GSVPoseOptimization::
Solve(RNBoolean update_translations, RNBoolean update_rotations, 
  RNBoolean update_laser_transformations, RNBoolean update_camera_transformations, 
  RNBoolean update_path_transformations, RNBoolean update_pixel_depths,
  RNBoolean include_scan_scan_correspondences, 
  RNBoolean include_image_image_correspondences, 
  RNBoolean include_scan_image_correspondences,
  RNScalar rigidity, int solver)
{
  // Determine which laser variables to optimize
  laser_nv = 0;
  for (int i = 0; i < LASER_NV; i++) laser_v[i] = -1;
  if (update_translations && update_laser_transformations) {
    laser_v[TX] = laser_nv++;
    laser_v[TY] = laser_nv++;
    laser_v[TZ] = laser_nv++;
  }

  // Determine which camera variables to optimize
  camera_nv = 0;
  for (int i = 0; i < CAMERA_NV; i++) camera_v[i] = -1;
  if (update_translations && update_camera_transformations) {
    camera_v[TX] = camera_nv++;
    camera_v[TY] = camera_nv++;
    camera_v[TZ] = camera_nv++;
  }

  // Determine which pixel variables to optimize
  pixel_nv = 0;
  for (int i = 0; i < PIXEL_NV; i++) pixel_v[i] = -1;
  if (update_pixel_depths) {
    pixel_v[T] = pixel_nv++;
  }

  // Determine which path variables to optimize
  vertex_nv = 0;
  for (int i = 0; i < VERTEX_NV; i++) vertex_v[i] = -1;
  if (update_translations && update_path_transformations) {
    vertex_v[TX] = vertex_nv++;
    vertex_v[TY] = vertex_nv++;
    vertex_v[TZ] = vertex_nv++;
  }
  if (update_rotations && update_path_transformations) {
    vertex_v[RX] = vertex_nv++;
    vertex_v[RY] = vertex_nv++;
    vertex_v[RZ] = vertex_nv++;
  }

  // Recompute list of image features ???
  pixels.Empty();
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    feature->pixel_index = -1;
  }
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    for (int j = 0; j < 2; j++) {
      GSVFeature *feature = correspondence->Feature(j);
      if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
      if (feature->pixel_index >= 0) continue;
      feature->pixel_index = pixels.NEntries();
      pixels.Insert(feature);
    }
  }

  // Get total number of variables
  int n = 0;
  n += laser_nv * lasers.NEntries();
  n += camera_nv * cameras.NEntries();
  n += vertex_nv * vertices.NEntries();
  n += pixel_nv * pixels.NEntries();
  if (n == 0) return 1;

  // Create system of equations
  RNPolynomialSystemOfEquations *equations = new RNPolynomialSystemOfEquations(n);
  if (!equations) {
    fprintf(stderr, "Unable to allocate system of equations\n");
    return 0;
  }
  
  printf("HEREA\n");

  // Insert equations 
  AddInertiaEquations(equations);
  printf("  HEREA1\n");
  AddRigidityEquations(equations, rigidity);
  printf("  HEREA2\n");
  AddCorrespondenceEquations(equations, include_scan_scan_correspondences, 
    include_image_image_correspondences, include_scan_image_correspondences);

  printf("HEREB\n");

  // Allocate x variables
  assert(n == equations->NVariables());
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = 0;

  // Initialize optimization variables
  if (!InitializeOptimizationVariables(x)) {
    fprintf(stderr, "Unable to initialize optimization variables\n");
    return 0;
  }

  // Optimize
  if (!ExecuteOptimization(equations, x, solver)) {
    fprintf(stderr, "Unable to extract optimization variables\n");
    return 0;
  }

  // Extract optimization variables
  if (!ExtractOptimizationVariables(x)) {
    fprintf(stderr, "Unable to extract optimization variables\n");
    return 0;
  }

  printf("HEREC\n");

  // Print errors
  printf("Final Errors:\n");
  PrintErrors(equations, x, rigidity, solver, TRUE);
  printf("\n");

  // Delete everything
  delete equations;
  delete [] x;

  // Return success
  return 1;
}



int GSVPoseOptimization::
ClearOptimizationVariables(void)
{
  // Clear all laser transformations
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    laser_data->translation = R3zero_vector;
  }

  // Clear all camera transformations
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    camera_data->translation = R3zero_vector;
  }

  // Clear all spline vertex transformations
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);
    vertex->translation = R3zero_vector;
    vertex->rotation = R3zero_vector;
  }

  // Clear all pixel depths
  for (int i = 0; i < pixels.NEntries(); i++) {
    const double default_image_t = 10;
    GSVFeature *feature = pixels.Kth(i);
    R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
    feature->scan_position = ray.Point(default_image_t);
    feature->image_t = 0;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
InitializeOptimizationVariables(double *x)
{
  // Initialize all laser variables
  if (laser_nv > 0) {
    for (int i = 0; i < lasers.NEntries(); i++) {
      for (int j = 0; j < LASER_NV; j++) {
        int v = LaserVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = LaserVariableValue(i, j);
      }
    }
  }

  // Initialize all camera variables
  if (camera_nv > 0) {
    for (int i = 0; i < cameras.NEntries(); i++) {
      for (int j = 0; j < CAMERA_NV; j++) {
        int v = CameraVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = CameraVariableValue(i, j);
      }
    }
  }

  // Initialize all vertex variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      for (int j = 0; j < VERTEX_NV; j++) {
        int v = VertexVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = VertexVariableValue(i, j);
      }
    }
  }

  // Initialize all pixel variables
  if (pixel_nv > 0) {
    for (int i = 0; i < pixels.NEntries(); i++) {
      for (int j = 0; j < PIXEL_NV; j++) {
        int v = PixelVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = PixelVariableValue(i, j);
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExtractOptimizationVariables(double *x)
{
  // Extract all laser variables
  if (laser_nv > 0) {
    printf("Lasers:\n");
    for (int i = 0; i < lasers.NEntries(); i++) {
      for (int j = 0; j < LASER_NV; j++) {
        SetLaserVariableValue(i, j, x);
        printf("%12.9f ", LaserVariableValue(i, j, x));
      }
      printf("\n");
    }
  }

  // Extract all camera variables
  if (camera_nv > 0) {
    printf("Cameras:\n");
    for (int i = 0; i < cameras.NEntries(); i++) {
      for (int j = 0; j < CAMERA_NV; j++) {
        SetCameraVariableValue(i, j, x);
        printf("%12.9f ", CameraVariableValue(i, j, x));
      }
      printf("\n");
    }
  }

  // Extract all vertex variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      for (int j = 0; j < VERTEX_NV; j++) {
        SetVertexVariableValue(i, j, x);
      }
    }
  }

  // Extract all pixel variables
  if (pixel_nv > 0) {
    for (int i = 0; i < pixels.NEntries(); i++) {
      for (int j = 0; j < PIXEL_NV; j++) {
        SetPixelVariableValue(i, j, x);
      }
    }
  }

  // Update pixel positions
  UpdatePixelPositionsFromCorrespondenceSolutions();

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExecuteOptimization(RNPolynomialSystemOfEquations *equations, double *x, int solver)
{
  // Solve system of equations
  if (!equations->Minimize(x, solver)) {
    fprintf(stderr, "Unable to minimize system of equations\n");
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Variable access utility functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
LaserVariableIndex(int laser_index, int variable_index)
{
  // Return index of laser variable
  if (laser_nv == 0) return -1;
  if (laser_v[variable_index] == -1) return -1;
  return laser_index*laser_nv + laser_v[variable_index];
}



int GSVPoseOptimization::
CameraVariableIndex(int camera_index, int variable_index)
{
  // Return index of camera variable
  if (camera_nv == 0) return -1;
  if (camera_v[variable_index] == -1) return -1;
  return lasers.NEntries()*laser_nv + 
         camera_index*camera_nv + camera_v[variable_index];
}



int GSVPoseOptimization::
VertexVariableIndex(int vertex_index, int variable_index)
{
  // Return index of spline variable
  if (vertex_nv == 0) return -1;
  if (vertex_v[variable_index] == -1) return -1;
  return lasers.NEntries()*laser_nv + 
         cameras.NEntries()*camera_nv + 
         vertex_index*vertex_nv + vertex_v[variable_index];
}



int GSVPoseOptimization::
PixelVariableIndex(int pixel_index, int variable_index)
{
  // Return index of pixel variable
  if (pixel_nv == 0) return -1;
  if (pixel_v[variable_index] == -1) return -1;
  return lasers.NEntries()*laser_nv + 
         cameras.NEntries()*camera_nv + 
         vertices.NEntries()*vertex_nv + 
         pixel_index*pixel_nv + pixel_v[variable_index];
}



////////////////////////////////////////



RNScalar GSVPoseOptimization::
LaserVariableValue(int laser_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (laser_v[variable_index] == -1)) 
    index = LaserVariableIndex(laser_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVLaser *laser = lasers.Kth(laser_index);
  LaserData *laser_data = (LaserData *) laser->Data();
  return laser_data->translation[variable_index];
}




RNScalar GSVPoseOptimization::
CameraVariableValue(int camera_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (camera_v[variable_index] == -1)) 
    index = CameraVariableIndex(camera_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVCamera *camera = cameras.Kth(camera_index);
  CameraData *camera_data = (CameraData *) camera->Data();
  return camera_data->translation[variable_index];
}




RNScalar GSVPoseOptimization::
VertexVariableValue(int vertex_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (vertex_v[variable_index] == -1)) 
    index = VertexVariableIndex(vertex_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVPathVertex *vertex = vertices.Kth(vertex_index);
  if (variable_index < 3) return vertex->translation[variable_index];
  else return vertex->rotation[variable_index-3];
}



RNScalar GSVPoseOptimization::
PixelVariableValue(int pixel_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (pixel_v[variable_index] == -1)) 
    index = PixelVariableIndex(pixel_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVFeature *pixel = pixels.Kth(pixel_index);
  if (variable_index == T) return pixel->image_t;
  else return 0.0;
}



////////////////////////////////////////



int GSVPoseOptimization::
SetLaserVariableValue(int laser_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (laser_v[variable_index] == -1) return 0; 
  int index = LaserVariableIndex(laser_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVLaser *laser = lasers.Kth(laser_index);
  LaserData *laser_data = (LaserData *) laser->Data();
  laser_data->translation[variable_index] = x[index];

  // Return success
  return 1;
}




int GSVPoseOptimization::
SetCameraVariableValue(int camera_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (camera_v[variable_index] == -1) return 0; 
  int index = CameraVariableIndex(camera_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVCamera *camera = cameras.Kth(camera_index);
  CameraData *camera_data = (CameraData *) camera->Data();
  camera_data->translation[variable_index] = x[index];

  // Return success
  return 1;
}




int GSVPoseOptimization::
SetVertexVariableValue(int vertex_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (vertex_v[variable_index] == -1) return 0; 
  int index = VertexVariableIndex(vertex_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVPathVertex *vertex = vertices.Kth(vertex_index);
  if (variable_index < 3) vertex->translation[variable_index] = x[index];
  else vertex->rotation[variable_index-3] = x[index];

  // Return success
  return 1;
}



int GSVPoseOptimization::
SetPixelVariableValue(int pixel_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (pixel_v[variable_index] == -1) return 0; 
  int index = PixelVariableIndex(pixel_index, variable_index);
  if (index == -1) return 0;

  // Return stored value
  GSVFeature *pixel = pixels.Kth(pixel_index);
  if (variable_index == T) pixel->image_t = x[index];

  // Return success
  return 1;
}



void GSVPoseOptimization::
PrintErrors(RNPolynomialSystemOfEquations *equations, double *x, double rigidity, int solver, int print_details, int print_values, int print_residuals)
{
  // Allocate memory
  int n = equations->NVariables();
  int m = equations->NPolynomials();
  double *y = new double [ m ];

  // Print parameters
  printf("Errors: ");
  printf("( %d %g | %d %d %d %d | %d %d %d %d | %d %d | %d %d %d %d %g )\n", 
    solver, rigidity, laser_nv, camera_nv, vertex_nv, pixel_nv, 
    lasers.NEntries(), cameras.NEntries(), vertices.NEntries(), pixels.NEntries(), 
    NFeatures(), NCorrespondences(), 
    equations->NVariables(), equations->NPolynomials(), 
    equations->NTerms(), equations->NPartialDerivatives(),
    equations->Degree());

  if (print_details) {
    // Print inertia error
    RNScalar inertia_error = 0;
    RNPolynomialSystemOfEquations inertia_equations(n);
    AddInertiaEquations(&inertia_equations);
    inertia_equations.Evaluate(x, y);
    int inertia_m = inertia_equations.NPolynomials();
    for (int i = 0; i < inertia_m; i++) inertia_error += y[i] * y[i];
    RNScalar inertia_rmsd = (inertia_m > 0) ? sqrt(inertia_error / inertia_m) : 0;
    printf("  %20s = %12.6f %12.6f %6d\n", "Inertia Error", inertia_error, inertia_rmsd, inertia_m);

    // Print rigidity error
    RNScalar rigidity_error = 0;
    RNPolynomialSystemOfEquations rigidity_equations(n);
    AddRigidityEquations(&rigidity_equations, rigidity);
    rigidity_equations.Evaluate(x, y);
    int rigidity_m = rigidity_equations.NPolynomials();
    for (int i = 0; i < rigidity_m; i++) rigidity_error += y[i] * y[i];
    RNScalar rigidity_rmsd = (rigidity_m > 0) ? sqrt(rigidity_error / rigidity_m) : 0;
    printf("  %20s = %12.6f %12.6f %6d\n", "Rigidity Error", rigidity_error, rigidity_rmsd, rigidity_m);

    // Print correspondence error
    RNScalar correspondence_error = 0;
    RNPolynomialSystemOfEquations correspondence_equations(n);
    AddCorrespondenceEquations(&correspondence_equations);   
    correspondence_equations.Evaluate(x, y);
    int correspondence_m = correspondence_equations.NPolynomials();
    for (int i = 0; i < correspondence_m; i++) correspondence_error += y[i] * y[i]; 
    RNScalar correspondence_rmsd = (correspondence_m > 0) ? sqrt(correspondence_error / correspondence_m) : 0;
    printf("  %20s = %12.6f %12.6f %6d\n", "Correspondence Error", correspondence_error, correspondence_rmsd, correspondence_m);
  }

  // Print overall error
  RNScalar overall_error = 0;
  equations->Evaluate(x, y);
  int overall_m = equations->NPolynomials();
  for (int i = 0; i < overall_m; i++) overall_error += y[i] * y[i];
  RNScalar overall_rmsd = (overall_m > 0) ? sqrt(overall_error / overall_m) : 0;
  printf("----------------------------------------------------------\n");
  printf("  %20s = %12.6f %12.6f %6d\n", "Overall Error", overall_error, overall_rmsd, overall_m);

  // Delete memory
  delete [] y;

  // Print values
  if (print_values) {
    printf("Values:\n");
    for (int i = 0; i < equations->NVariables(); i++) {
      printf("  %9.6f ", x[i]);
      if ((i > 0) && ((i % 6) == 0)) printf("\n");
    }
    printf("\n");
  }

  // Print residuals
  if (print_residuals) {
    int m = equations->NPolynomials();
    double *y = new double [ m ];
    equations->Evaluate(x, y);
    printf("Residuals:\n");
    for (int i = 0; i < m; i++) 
      printf("  %9.6f\n", y[i]);
    printf("\n");
    delete [] y;
  }
}

