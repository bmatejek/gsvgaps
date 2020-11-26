// Include file for pose optimization 



////////////////////////////////////////////////////////////////////////
// Dependencies
////////////////////////////////////////////////////////////////////////

#include "RNMath/RNMath.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

struct GSVDescriptor {
public:
  GSVDescriptor(int descriptor_type = 0, RNScalar *values = NULL, int nvalues = 0);
  GSVDescriptor(const GSVDescriptor& descriptor);
  ~GSVDescriptor(void);
  RNScalar SquaredDistance(const GSVDescriptor& descriptor) const;
public:
  int descriptor_type;
  RNScalar *values; 
  int nvalues;
};

struct GSVFeature {
public:
  GSVFeature(void);
  GSVFeature(int feature_type, GSVScanline *scanline, int scan_point_index, 
    const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, 
    RNLength scan_scale, RNScalar score);
  GSVFeature(int feature_type, GSVImage *image, 
    const R2Point& image_position, const R2Vector& image_direction, 
    RNLength image_scale, RNScalar score);
  GSVFeature(int feature_type, GSVImage *image, GSVScanline *scanline, int scan_point_index, 
    const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, RNLength scan_scale,
    const R2Point& image_position, const R2Vector& image_direction, RNLength image_scale, 
    RNScalar image_t, RNScalar score);
  void Draw(void) const;
public:
  int feature_type;
  GSVImage *image;
  GSVScanline *scanline;
  int scan_point_index;
  R3Point scan_position;
  R3Vector scan_direction;
  R3Vector scan_normal;
  RNScalar scan_scale;
  R2Point image_position;
  R2Vector image_direction;
  RNScalar image_scale;
  RNScalar image_t;
  RNScalar score;
  GSVDescriptor descriptor;
  int pixel_index;
  int index;
};

struct GSVFeatureSet {
public:
  GSVFeatureSet(GSVScene *scene);
  GSVFeatureSet(const GSVFeatureSet& set);
  int NFeatures(void) const;
  GSVFeature *Feature(int k) const;
  GSVFeatureSet& operator=(const GSVFeatureSet& set);
  void InsertFeature(GSVFeature *feature);
  void RemoveFeature(GSVFeature *feature);
  void Empty(void);
  void Draw(void) const;
public:
  GSVScene *scene;
  RNArray<GSVFeature *> features;
};

struct GSVFeatureCorrespondence {
public:
  GSVFeatureCorrespondence(void);
  GSVFeatureCorrespondence(GSVFeature *feature0, GSVFeature *feature1, RNScalar score = 1);
  GSVFeature *Feature(int k) const;
  void Draw(void) const;
public:
  GSVFeature *features[2];
  RNScalar score;
  int index;
};

struct GSVFeatureCorrespondenceSet {
public:
  GSVFeatureCorrespondenceSet(GSVScene *scene);
  GSVFeatureCorrespondenceSet(const GSVFeatureCorrespondenceSet& set);
  int NCorrespondences(void) const;
  GSVFeatureCorrespondence *Correspondence(int k) const;
  GSVFeatureCorrespondenceSet& operator=(const GSVFeatureCorrespondenceSet& set);
  void InsertCorrespondence(GSVFeatureCorrespondence *correspondence);
  void RemoveCorrespondence(GSVFeatureCorrespondence *correspondence);
  void Empty(void);
  void Sort(void);
  void Draw(void) const;
public:
  GSVScene *scene;
  RNArray<GSVFeatureCorrespondence *> correspondences;
  RNScalar score;
};

class GSVPathVertex {
public:
  GSVPathVertex(void);
public:
  struct GSVPath *path;
  GSVPose pose;
  R3Vector translation;
  R3Vector rotation;
  RNScalar parameter;
  RNScalar timestamp;
  int index;
};

struct GSVPath {
public:
  GSVPath(GSVSegment *segment = NULL, RNLength max_vertex_spacing = 2);
  ~GSVPath(void);
  GSVPose TransformedPose(RNScalar u) const;
  GSVPose Pose(RNScalar u) const;
  R3Vector Translation(RNScalar u) const;
  R3Vector Rotation(RNScalar u) const;
  RNScalar Parameter(RNScalar timestamp) const;
  RNScalar Timestamp(RNScalar u) const;
  RNInterval ParameterRange(void) const;
  RNInterval TimestampRange(void) const;
  void Draw(void) const;
public:
  friend struct GSVPoseOptimization;
  int NVertices(void) const;
  GSVPathVertex *Vertex(int vertex_index) const;
  const GSVPose& VertexPose(int vertex_index) const;
  const R3Point& VertexPosition(int vertex_index) const;
  const R3Vector& VertexTranslation(int vertex_index) const;
  const R3Vector& VertexRotation(int vertex_index) const;
  RNScalar VertexParameter(int vertex_index) const;
  RNScalar VertexTimestamp(int vertex_index) const;
  GSVPose VertexTransformedPose(int vertex_index) const;
  void SetVertexTranslation(int vertex_index, const R3Vector& rotation);
  void SetVertexRotation(int vertex_index, const R3Vector& rotation);
  int FindVertexIndexBeforeTimestamp(RNScalar timestamp) const;
  int FindVertexIndexBeforeTimestamp(RNScalar timestamp, int imin, int imax) const;
public:
  GSVSegment *segment;
  R3CatmullRomSpline *spline;
  RNArray<class GSVPathVertex *> vertices;
  int index;
};

struct GSVPoseOptimization {
public:
  GSVPoseOptimization(GSVScene *scene);
  RNScalar Score(void);
  RNScalar RMSD(void);
  int NFeatures(void) const;
  GSVFeature *Feature(int k) const;
  int NCorrespondences(void) const;
  GSVFeatureCorrespondence *Correspondence(int k) const;
  R3Affine OptimizedTransformation(const GSVScanline *scanline) const;
  R3Affine OptimizedTransformation(const GSVImage *image) const;
  R3Point FeaturePosition(const GSVFeature *feature, RNBoolean optimized = FALSE) const;
  R3Vector FeatureDirection(const GSVFeature *feature, RNBoolean optimized = FALSE) const;
  R3Vector FeatureNormal(const GSVFeature *feature, RNBoolean optimized = FALSE) const;
  RNRgb FeatureColor(GSVFeature *feature) const;
  void InsertFeature(GSVFeature *feature);
  void RemoveFeature(GSVFeature *feature);
  void InsertCorrespondence(GSVFeatureCorrespondence *correspondence);
  void RemoveCorrespondence(GSVFeatureCorrespondence *correspondence);
  int CreateImageCornerFeatures(void);
  int CreateImageSiftFeatures(void);
  int CreateImageLineFeatures(void);
  int CreateScanPoleFeatures(void);
  int CreateScanCurbFeatures(void);
  int CreateScanEdgeFeatures(void);
  int CreateScanPlaneFeatures(void);
  int CreateCorrespondences(void);
  int ReadFile(const char *filename);
  int WriteFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadAsciiFile(const char *filename);
  int WriteAsciiFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadAscii(FILE *fp);
  int WriteAscii(FILE *fp, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadBinaryFile(const char *filename);
  int WriteBinaryFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadBinary(FILE *fp);
  int WriteBinary(FILE *fp, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadSiftFile(GSVImage *image, const char *filename);
  int ReadSiftFile(const char *filename);
  int ReadLineFile(GSVImage *image, const char *filename);
  int ReadLineFile(const char *filename);
  int ReadTransformationsFile(const char *filename);
  int WriteTransformationsFile(const char *filename) const;
  int ApplyOptimizedTransformationsToScene(void);
  int Solve(
    RNBoolean update_translations = TRUE, RNBoolean update_rotations = TRUE, 
    RNBoolean update_laser_transformations = TRUE, RNBoolean update_camera_transformations = TRUE, 
    RNBoolean update_path_transformations = TRUE, RNBoolean update_pixel_depths = TRUE,
    RNBoolean include_scan_scan_correspondences = TRUE, 
    RNBoolean include_image_image_correspondences = TRUE, 
    RNBoolean include_scan_image_correspondences = TRUE,
    RNScalar rigidity = 0.5, int solver = 0);
public:
  RNScalar CorrespondenceScore(const GSVFeature *feature1, const GSVFeature *feature2);
  RNBoolean IsMatchFeasible(const GSVFeature *feature1, const GSVFeature *feature2, const GSVFeatureCorrespondenceSet *match = NULL);
  int CreateAllFeasibleCorrespondences(void);
  int CreateCoplanarCorrespondences(void);
  int CreateImageToScanCorrespondences(void);
  int CreateMutuallyClosestCorrespondences(void);
  int ExhaustiveSearch(const GSVFeatureCorrespondenceSet& pairwise_matches, 
    int max_correspondences, GSVFeatureCorrespondenceSet& best_match);
  int ExhaustiveSearch(void);
  int RANSAC(void);
  int ICP(void);
  int UpdatePixelPositionsFromCorrespondenceSolutions(void);
  int UpdatePixelPositionsFromRayIntersections(void);
  int ClearPixelPositions(void);
public:
  GSVScene *scene;
  GSVFeatureSet features;
  GSVFeatureCorrespondenceSet correspondences;
  RNArray<GSVLaser *> lasers;
  RNArray<GSVCamera *> cameras;
  RNArray<GSVSegment *> segments;
  RNArray<GSVScanline *> scanlines;
  RNArray<GSVImage *> images;
  RNArray<GSVPathVertex *> vertices;
  RNScalar score;

  // Feature parameters
  RNScalar sift_image_scale;

  // Correspondence parameters
  RNBoolean create_point_point_correspondences;
  RNBoolean create_plane_plane_correspondences;
  RNBoolean create_pixel_point_correspondences;
  RNBoolean create_pixel_plane_correspondences;
  RNBoolean create_pixel_pixel_correspondences;
  RNScalar min_intracorrespondence_path_parameter_difference;
  RNScalar max_intracorrespondence_euclidean_distance;
  RNScalar max_intracorrespondence_spin_image_descriptor_distance;
  RNScalar max_intracorrespondence_sift_descriptor_distance;
  RNScalar max_intracorrespondence_line_descriptor_distance;
  RNScalar max_intracorrespondence_direction_angle;
  RNScalar max_intracorrespondence_normal_angle;
  RNScalar min_intrasegment_path_parameter_difference;
  RNScalar max_intrasegment_distance_ratio;
  RNScalar icp_max_distance_start;
  RNScalar icp_max_distance_end;

  // Optimization variable indicies
  static const int LASER_NV = 3;
  static const int CAMERA_NV = 3;
  static const int VERTEX_NV = 6;
  static const int TX = 0;
  static const int TY = 1;
  static const int TZ = 2;
  static const int RX = 3;
  static const int RY = 4;
  static const int RZ = 5;
  static const int PIXEL_NV = 1;
  static const int T = 0;

  // Optimization variables
  RNArray<GSVFeature *> pixels;
  int laser_v[LASER_NV];
  int laser_nv;
  int camera_v[CAMERA_NV];
  int camera_nv;
  int vertex_v[VERTEX_NV];
  int vertex_nv;
  int pixel_v[PIXEL_NV];
  int pixel_nv;

  // Optimization parameters
  RNScalar laser_inertia_weights[LASER_NV];
  RNScalar camera_inertia_weights[CAMERA_NV];
  RNScalar vertex_inertia_weights[VERTEX_NV];
  RNScalar pixel_inertia_weights[PIXEL_NV];
  int intra_segment_rigidity_radius;
  RNScalar intra_segment_rigidity_sigma;
  RNScalar intra_segment_rigidity_sigma_squared;
  RNScalar intra_segment_rigidity_weight;
  RNScalar scan_point_scan_point_correspondence_weight;
  RNScalar scan_point_scan_plane_correspondence_weight;
  RNScalar scan_point_image_point_correspondence_weight;
  RNScalar image_point_scan_plane_correspondence_weight;
  RNScalar image_point_image_point_correspondence_weight;
  RNScalar max_vertex_spacing;

  // OPTIMIZATION STUFF

  // Optimization functions
  int ClearOptimizationVariables(void);
  int InitializeOptimizationVariables(double *x);
  int ExtractOptimizationVariables(double *x);
  int ExecuteOptimization(RNPolynomialSystemOfEquations *system, double *x, int solver);

  // Inertial equations
  void AddInertiaEquations(RNPolynomialSystemOfEquations *system);

  // Rigidity equations
  void AddRigidityEquations(RNPolynomialSystemOfEquations *system, RNScalar rigidity); 
  void AddPathRigidityEquations(RNPolynomialSystemOfEquations *system, 
    int spline_index0, int spline_index1, const R3Point& v0, const R3Point& v1, RNScalar w);

  // Correspondence equations
  void AddCorrespondenceEquations(RNPolynomialSystemOfEquations *system,
    RNBoolean include_scan_scan_correspondences = TRUE, 
    RNBoolean include_image_image_correspondences = TRUE, 
    RNBoolean include_scan_image_correspondences = TRUE);
  void AddCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVFeatureCorrespondence *correspondence,
    RNBoolean include_scan_scan_correspondences = TRUE, 
    RNBoolean include_image_image_correspondences = TRUE, 
    RNBoolean include_scan_image_correspondences = TRUE);
  void AddPointPointCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVScanline *scanline0, const R3Point& point0, 
    GSVScanline *scanline1, const R3Point& point1, 
    RNScalar weight);
  void AddPointRayCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVScanline *scanline0, const R3Point& point0, 
    GSVImage *image1, const R3Ray& ray1, 
    RNScalar weight);
  void AddPointPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVScanline *scanline0, const R3Point& point0, 
    GSVScanline *scanline1, const R3Point& point1, const R3Vector& normal1, 
    RNScalar weight);
  void AddPointPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVScanline *scanline0, const R3Point& point0, 
    GSVImage *image1, const R3Point& point1, const R3Vector& normal1, 
    RNScalar weight);
  void AddPixelPointCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVImage *image0, const R2Point& point0, int pixel_index0, 
    GSVScanline *scanline1, const R3Point& point1, 
    RNScalar weight);
  void AddPixelRayCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVImage *image0, const R2Point& point0, int pixel_index0, 
    GSVImage *image1, const R3Ray& ray1, 
    RNScalar weight);
  void AddPixelPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVImage *image0, const R2Point& point0, int pixel_index0, 
    GSVScanline *scanline1, const R3Point& point1, const R3Vector& normal1, 
    RNScalar weight);
  void AddPixelPlaneCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVImage *image0, const R2Point& point0, int pixel_index0, 
    GSVImage *image1, const R3Point& point1, const R3Vector& normal1, 
    RNScalar weight);
  void AddPixelPixelCorrespondenceEquations(RNPolynomialSystemOfEquations *system, 
    GSVImage *image0, const R2Point& point0, int pixel_index0, 
    GSVImage *image1, const R2Point& pixel1, int pixel_index1, 
    RNScalar weight);

  // Mid-level system of equation utility functions
  int ComputeTransformedPointCoordinates(GSVScanline *scanline, const R3Point& point,
    RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz);
  int ComputeTransformedPointCoordinates(GSVImage *image, const R3Point& point,
    RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz);
  int ComputeTransformedPointCoordinates(GSVImage *image, const R2Point& point, int pixel_index,
    RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz);
  int ComputeTransformedVectorCoordinates(GSVScanline *scanline, const R3Vector& vector,
    RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz);
  int ComputeTransformedVectorCoordinates(GSVImage *image, const R3Vector& vector,
    RNPolynomial& px, RNPolynomial& py, RNPolynomial& pz);

  // Low-level system of equation utility functions
  void AddLaserVariableToEquation(RNPolynomial *equation, int laser_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddCameraVariableToEquation(RNPolynomial *equation, int camera_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPathVertexVariableToEquation(RNPolynomial *equation, int vertex_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPathVertexVariableToEquation(RNPolynomial *equation, GSVPath *path, RNScalar path_parameter, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddScanlineVariableToEquation(RNPolynomial *equation, int scanline_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddImageVariableToEquation(RNPolynomial *equation, int image_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPixelVariableToEquation(RNPolynomial *equation, int pixel_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);

  // Variable index utilities
  int LaserVariableIndex(int pixel_index, int variable_index);
  int CameraVariableIndex(int pixel_index, int variable_index);
  int VertexVariableIndex(int pixel_index, int variable_index);
  int PixelVariableIndex(int pixel_index, int variable_index);

  // Variable value utilities
  RNScalar LaserVariableValue(int laser_index, int variable_index, RNScalar *x = NULL);
  RNScalar CameraVariableValue(int camera_index, int variable_index, RNScalar *x = NULL);
  RNScalar VertexVariableValue(int vertex_index, int variable_index, RNScalar *x = NULL);
  RNScalar PixelVariableValue(int pixel_index, int variable_index, RNScalar *x = NULL);

  // Variable assignment utilities
  int SetLaserVariableValue(int laser_index, int variable_index, RNScalar *x);
  int SetCameraVariableValue(int camera_index, int variable_index, RNScalar *x);
  int SetVertexVariableValue(int vertex_index, int variable_index, RNScalar *x);
  int SetPixelVariableValue(int pixel_index, int variable_index, RNScalar *x);

  // Error reporting
  void PrintErrors(RNPolynomialSystemOfEquations *equations, double *x, double rigidity, int solver, 
    int print_details = 0, int print_values = 0, int print_residuals = 0);
};



////////////////////////////////////////////////////////////////////////
// Other type defs???
////////////////////////////////////////////////////////////////////////

struct LaserData {
  GSVLaser *laser;
  R3Vector translation;
  int index;
};

struct CameraData {
  GSVCamera *camera;
  R3Vector translation;
  int index;
};

struct SegmentData {
  GSVSegment *segment;
  GSVPath *path;
  int index;
};

struct ScanData {
  GSVScan *scan;
  GSVMesh *mesh;
  R2Grid da_scanline_grid;
  R2Grid da_point_index_grid;
  R2Grid da_position_x_grid;
  R2Grid da_position_y_grid;
  R2Grid da_position_z_grid;
  R2Grid da_plane_segment_id_grid;
  RNScalar timestamp;
  int index;
};

struct ImageData {
  GSVImage *image;
  RNScalar path_parameter;
  int index;
};

struct ScanlineData {
  GSVScanline *scanline;
  RNScalar path_parameter;
  int index;
};



////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////

enum {
  GSV_NULL_DESCRIPTOR_TYPE,
  GSV_SPIN_IMAGE_DESCRIPTOR_TYPE,
  GSV_SIFT_DESCRIPTOR_TYPE,
  GSV_SHAPE_CONTEXT_DESCRIPTOR_TYPE,
  GSV_LINE_DESCRIPTOR_TYPE,
  GSV_NUM_DESCRIPTOR_TYPES
};

enum {
  GSV_NULL_FEATURE_TYPE,
  GSV_SCAN_POINT_FEATURE_TYPE,
  GSV_SCAN_PLANE_FEATURE_TYPE,
  GSV_SCAN_LINE_FEATURE_TYPE,
  GSV_IMAGE_POINT_FEATURE_TYPE,
  GSV_IMAGE_LINE_FEATURE_TYPE,
  GSV_NUM_FEATURE_TYPES
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline int GSVFeatureSet::
NFeatures(void) const
{
  // Return number of features
  return features.NEntries();
}



inline GSVFeature *GSVFeatureSet::
Feature(int k) const
{
  // Return kth feature
  return features.Kth(k);
}



inline GSVFeature *GSVFeatureCorrespondence::
Feature(int k) const
{
  // Return kth feature
  assert((k == 0) || (k == 1));
  return features[k];
}



inline int GSVFeatureCorrespondenceSet::
NCorrespondences(void) const
{
  // Return number of correspondences
  return correspondences.NEntries();
}



inline GSVFeatureCorrespondence *GSVFeatureCorrespondenceSet::
Correspondence(int k) const
{
  // Return kth correspondence
  return correspondences.Kth(k);
}



inline int GSVPoseOptimization::
NFeatures(void) const
{
  // Return number of features
  return features.NFeatures();
}



inline GSVFeature *GSVPoseOptimization::
Feature(int k) const
{
  // Return kth feature
  return features.Feature(k);
}



inline int GSVPoseOptimization::
NCorrespondences(void) const
{
  // Return number of correspodnences
  return correspondences.NCorrespondences();
}



inline GSVFeatureCorrespondence *GSVPoseOptimization::
Correspondence(int k) const
{
  // Return kth correspodnence
  return correspondences.Correspondence(k);
}



