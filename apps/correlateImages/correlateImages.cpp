////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include "correlateImages.h"
#include "findBBox.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Program Arguments
////////////////////////////////////////////////////////////////////////

// global macros
#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

// command line variables
static int print_verbose = 0;
static int run = -1;
static int segment = -1;
static double nearHitRange = 2;
static double gsvToVOCRatio = 1;
static char *input_scene_name = NULL;
static char *input_gsv_correspondences_name = NULL;
static char *output_data_name = NULL;
vector<char *> voc_files = vector<char *>();

////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////

static int indexFromIndicies(int tapestry_index, int panorama_index) {
	return 9 * tapestry_index + panorama_index;
}

static int indexFromImage(GSVImage *image) {
	return indexFromIndicies(image->TapestryIndex(), image->PanoramaIndex());
}

static int reverseTapestryLookup(int index) {
	return index / 9;
}

static int reversePanoramaLookup(int index) {
	return index % 9;
}

// convert run from normal indices (1, 2, 3) to gsv indices (1, 10, 11 ... 19, 2, ...)
static int convertRun(int run) {
  int run_number = run;
  if (run >= 2 && run <= 9) run_number += 10;
  else if (run >= 10 && run <= 19) run_number -= 8;
  return run_number;
}

// convert GSV Points to VOC points
static R2Point ConvertGSVToVOC(R2Point gsv) {
  return R2Point(gsv.X() / gsvToVOCRatio, (GSV_HEIGHT - gsv.Y()) / gsvToVOCRatio);
}	

// return true if the point is within bbox and false otherwise
static bool InsideBBox(R2Point point, R2Box bbox) {
  return bbox.XMin() < point.X() && bbox.YMin() < point.Y() && point.X() < bbox.XMax() && point.Y() < bbox.YMax();
}

// return true if the point is near bbox and false otherwise
static bool NearBBox(R2Point point, R2Box bbox) {
  R2Box inflated = R2Box(bbox);
  inflated.Inflate(nearHitRange);
  return InsideBBox(point, inflated);
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

// read in the GSVScene
static GSVScene *ReadScene(const char *filename) {
  // start statistics
  RNTime start_time;
  start_time.Read();
  
  // allocate google scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }
  
  // read scene
  if (!scene->ReadFile(filename, 0)) {
    delete scene;
    return NULL;
  }
  
  // print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
  // return scene
  return scene;
}

// read the correspondences between lidar and image features
static int ReadCorrespondences(GSVPoseOptimization *optimization, const char *filename) {
  // start statistics
  RNTime start_time;
  start_time.Read();
  
  // read correspondences
  if (!optimization->ReadFile(filename)) return 0;
  
  // print statistics
  if (print_verbose) {
    printf("Read correspondenecs from %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
  // return OK status
  return 1;
}

// read VOC features into a list
static int ReadVOCFeatures(GSVPoseOptimization *optimization, VOCFeatureList *vocFeatureList, const char *filename) {
  // start statistics
  RNTime start_time;
  start_time.Read();
  
  // THRESHOLD HARDCODED IN
  double threshold = 0.5;
  
  // read correspondences
  if (!vocFeatureList->ReadFile(filename, optimization, threshold)) return 0;
  
  // print statistics
  if (print_verbose) {
    printf("Read VOC features from %s ...\n", filename);
    printf("\tTime = %.2f seconds\n", start_time.Elapsed());
  }
  
  // return OK status
  return 1;
}

////////////////////////////////////////////////////////////////////////
// CorrelateImage Constructor
////////////////////////////////////////////////////////////////////////

CorrelateImage ::
CorrelateImage(int run, int segment, char *input_scene_name, char *input_gsv_correspondences_name) :
	run(run),
	segment(segment),
	optimization(NULL),
	images(vector<Image>()),
	numberOfImages(0)
{
	// read input scene
	GSVScene *scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);
	
	// allocate optimization data structure
	optimization = new GSVPoseOptimization(scene);
	if (!optimization) { delete scene; exit(-1); }
	
	// read correspondences
	if (!ReadCorrespondences(optimization, input_gsv_correspondences_name)) { delete scene; delete optimization; exit(-1); }

	run = convertRun(run);
	int numImages = scene->Run(run)->Segment(segment)->NImages();
	int numPanoramas = scene->Run(run)->Segment(segment)->NPanoramas();
	// initialize all of the VOC vectors in the image array
	for (int i = 0; i < numPanoramas; i++) {
		GSVPanorama *panorama = scene->Run(run)->Segment(segment)->Panorama(i);
		int numPanoramaImages = panorama->NImages();
		for (int j = 0; j < numPanoramaImages; j++) {
			Image addImage;
			addImage.vocFeatures = vector<VOCFeature *>();
			addImage.gsvFeatures = vector<GSVFeature *>();
			addImage.lidarImagePairs = vector<LidarImagePair *>();
			addImage.image = panorama->Image(j);
			images.push_back(addImage);
		}
	}
	
	// add gsv image features to images array
	for (int i = 0; i < optimization->NFeatures(); i++) {
		GSVFeature *feature = optimization->Feature(i);
		// make sure feature is image feature
		if (feature->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
			int index = indexFromImage(feature->image);
			images[index].gsvFeatures.push_back(feature);
		}
	}

	// add lidar image pairs to images array
	for (int i = 0; i < optimization->NCorrespondences(); i++) {
		GSVFeature *imageFeature = optimization->Correspondence(i)->features[0];
		GSVFeature *lidarFeature = optimization->Correspondence(i)->features[1];
		// only consider unique lidar points
		if (imageFeature->feature_type != lidarFeature->feature_type) {
			// switch lidar and image features
			if (imageFeature->feature_type == GSV_SCAN_POINT_FEATURE_TYPE && lidarFeature->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
				GSVFeature *temp = lidarFeature;
				lidarFeature = imageFeature;
				imageFeature = temp;
			}
			// see if the lidar feature already exists
			struct GSVLidarFeature *gsvLidarFeature = NULL;
			for (unsigned int j = 0; j  < gsvLidarFeatures.size(); j++) {
				if (gsvLidarFeatures[j]->lidar == lidarFeature) {
					gsvLidarFeature = gsvLidarFeatures[j];
				}
			}
			// if not, create the GSVLidarFeature
			if (!gsvLidarFeature) {
				gsvLidarFeature = (GSVLidarFeature *) malloc(sizeof(GSVLidarFeature));
				gsvLidarFeature->lidar = lidarFeature;
				gsvLidarFeature->matched = false;
				gsvLidarFeature->nearlyMatched = false;
				gsvLidarFeatures.push_back(gsvLidarFeature);
			}
			// create the gsv lidar image pair
			struct LidarImagePair *lidarImagePair = (LidarImagePair *) malloc(sizeof(LidarImagePair));
			lidarImagePair->gsvLidarFeature = gsvLidarFeature;
			lidarImagePair->gsvImageFeature = imageFeature;
			int index = indexFromImage(imageFeature->image);
			images[index].lidarImagePairs.push_back(lidarImagePair);
		}
	}
	numberOfImages = numImages;
}

////////////////////////////////////////////////////////////////////////
// CorrelateImage Accessors/Manipulators
////////////////////////////////////////////////////////////////////////

Image CorrelateImage ::
Images(int i) {
	return images[i];
}

GSVPoseOptimization *CorrelateImage ::
Optimization(void) {
	return optimization;
}

int CorrelateImage ::
AddVOCFeatures(char *input_voc_features_name) {

	// allocate VOCFeatureList data structure
	vocFeatures = new VOCFeatureList();
	if (!vocFeatures) { return 0; }
		
	// read VOC features
	if (!ReadVOCFeatures(optimization, vocFeatures, input_voc_features_name)) { return 0; }
	
	// clear all VOC vectors in the image array
	for (int i = 0; i < numberOfImages; i++) {
		images[i].vocFeatures.clear();
	}
	
	// reset all of the lidar matches
	for (unsigned int i = 0; i < gsvLidarFeatures.size(); i++) {
		gsvLidarFeatures[i]->matched = false;
		gsvLidarFeatures[i]->nearlyMatched = false;
	}
	
	// add voc features to images array
	for (int i = 0; i < vocFeatures->NFeatures(); i++) {
		VOCFeature *feature = (VOCFeature *) vocFeatures->Feature(i);
		int index = indexFromImage((GSVImage *) feature->Image());
		images[index].vocFeatures.push_back(feature);
	}
	
	// return OK status
	return 1;
}

const int CorrelateImage ::
NImages(void) const {
	return numberOfImages;
}

double CorrelateImage ::
Recall(void) {
  int totalGSVFeatures = 0;
  int totalHits = 0;
  int totalNearHits = 0;
  int totalLidarHits = 0;
  int totalLidarNearHits = 0;
  int totalLidarFeatures = (int) gsvLidarFeatures.size();
  
  for (int i = 0; i < numberOfImages; i++) {
    // only allow VOC bboxes to be used once
    bool vocUsed[images[i].vocFeatures.size()];
    bool vocNearlyUsed[images[i].vocFeatures.size()];
    for (unsigned int k = 0; k < images[i].vocFeatures.size(); k++) {
      vocUsed[k] = false;
      vocNearlyUsed[k] = false;
    }
    for (unsigned int j = 0; j < images[i].gsvFeatures.size(); j++) {
      totalGSVFeatures++;
      GSVFeature *gsvFeature = images[i].gsvFeatures[j];
      R2Point gsvPoint = ConvertGSVToVOC(gsvFeature->image_position);
      bool hit = false;
      bool nearHit = false;			
      for (unsigned int k = 0; k < images[i].vocFeatures.size(); k++) {
	VOCFeature *vocFeature = images[i].vocFeatures[k];
	// if the point is in the box, there is a hit
	if (!hit && !vocUsed[k] && InsideBBox(gsvPoint, vocFeature->BBox())) {
	  totalHits++;
	  hit = true;
	  vocUsed[k] = true;
	  // determine if there is a lidar match
	  for (unsigned int w = 0; w < images[i].lidarImagePairs.size(); w++) {
	    if (images[i].lidarImagePairs[w]->gsvImageFeature == gsvFeature) {
	      images[i].lidarImagePairs[w]->gsvLidarFeature->matched = true;
	    }
	  }
	}
	// if the point is near the box, there is a near hit
	if (!nearHit && !vocNearlyUsed[k] && NearBBox(gsvPoint, vocFeature->BBox())) {
	  totalNearHits++;
	  nearHit = true;
	  vocNearlyUsed[k] = true;
	  // determine if there is a lidar match
	  for (unsigned int w = 0; w < images[i].lidarImagePairs.size(); w++) {
	    if (images[i].lidarImagePairs[w]->gsvImageFeature == gsvFeature) {
	      images[i].lidarImagePairs[w]->gsvLidarFeature->nearlyMatched = true;
	    }
	  }
	}
      }
    }
  }
  
  for (unsigned int i = 0; i < gsvLidarFeatures.size(); i++) {
    if (gsvLidarFeatures[i]->matched)
      totalLidarHits++;
    if (gsvLidarFeatures[i]->nearlyMatched)
      totalLidarNearHits++;
  }
  
  if (print_verbose) {
    printf("\n");
    printf("Percent Image Recall:\n");
    printf("\tTrue Hits: %d/%d = %.3f%%\n", totalHits, totalGSVFeatures, 100.0 * ((double) totalHits) / ((double) totalGSVFeatures));
    printf("\tNear Hits: %d/%d = %.3f%%\n", totalNearHits, totalGSVFeatures, 100.0 * ((double) totalNearHits) / ((double) totalGSVFeatures));
    printf("\n");
    printf("Percent Lidar Recall:\n");
    printf("\tTrue Hits: %d/%d = %.3f%%\n", totalLidarHits, totalLidarFeatures, 100.0 * ((double) totalLidarHits) / ((double) totalLidarFeatures));
    printf("\tNear Hits: %d/%d = %.3f%%\n", totalLidarNearHits, totalLidarFeatures, 100.0 * ((double) totalLidarNearHits) / ((double) totalLidarFeatures));
  }
  
  return 100.0 * ((double) totalNearHits) / ((double) totalGSVFeatures);
}

double CorrelateImage ::
Precision(void) {
  int totalVOCFeatures = 0;
  int totalVOCHits = 0;
  int totalVOCNearHits = 0;
  
  for (int i = 0; i < numberOfImages; i++) {
    // make sure gsv points are only used once
    bool gsvUsed[images[i].gsvFeatures.size()];
    bool gsvNearlyUsed[images[i].gsvFeatures.size()];
    for (unsigned int k = 0; k < images[i].gsvFeatures.size(); k++) {
      gsvUsed[k] = false;
      gsvNearlyUsed[k] = false;
    }
    for (unsigned int j = 0; j < images[i].vocFeatures.size(); j++) {
      totalVOCFeatures++;
      bool hit = false;
      bool nearHit = false;
      VOCFeature *vocFeature = images[i].vocFeatures[j];
      for (unsigned int k = 0; k < images[i].gsvFeatures.size(); k++) {
	GSVFeature *gsvFeature = images[i].gsvFeatures[k];
	R2Point gsvPoint = ConvertGSVToVOC(gsvFeature->image_position);
	// if the point is inside the box, there is a hit
	if (!hit && !gsvUsed[k] && InsideBBox(gsvPoint, vocFeature->BBox())) {
	  totalVOCHits++;
	  hit = true;
	  gsvUsed[k] = true;
	}
	// if the point is near the box, there is a near hit
	if (!nearHit && !gsvNearlyUsed[k] && NearBBox(gsvPoint, vocFeature->BBox())) {
	  totalVOCNearHits++;
	  nearHit = true;
	  gsvNearlyUsed[k] = true;
	}
      }
    }
  }
  
  if (print_verbose) {
    printf("\n");
    printf("Percent VOC Precision:\n");
    printf("\tTrue Hits: %d/%d = %.3f%%\n", totalVOCHits, totalVOCFeatures, 100.0 * ((double) totalVOCHits) / ((double) totalVOCFeatures));
    printf("\tNear Hits: %d/%d = %.3f%%\n", totalVOCNearHits, totalVOCFeatures, 100.0 * ((double) totalVOCNearHits) / ((double) totalVOCFeatures));
  }
  return 100.0 * ((double) totalVOCNearHits) / ((double) totalVOCFeatures);
}

////////////////////////////////////////////////////////////////////////
// Output Functions
////////////////////////////////////////////////////////////////////////

void CorrelateImage ::
OutputData(void) {
  ofstream output;
  output.open(output_data_name, ios::out | ios::app);
  
  // return FAILURE
  if (!output) { fprintf(stderr, "Failed to write output file: %s\n", output_data_name); return; }
  
  for (unsigned int i = 0; i < images.size(); i++) {
    // only print if the image has a gsv feature and/or a voc feature
    if (images[i].gsvFeatures.size() || images[i].vocFeatures.size()) {
      int run_number = run;
      int segment_number = segment;
      int tapestry_index = reverseTapestryLookup(i);
      int panorama_index = reversePanoramaLookup(i);
      
      // output image information
      output << run_number; output << " ";
      output << segment_number; output << " ";
      output << tapestry_index; output << " ";
      output << panorama_index;
      
      // print out (x,y) coordinates of GSV features
      for (unsigned int j = 0; j < images[i].gsvFeatures.size(); j++) {
	R2Point gsvPoint = ConvertGSVToVOC(images[i].gsvFeatures[j]->image_position);
	output << " (";
	output << gsvPoint.X(); output << ",";
	output << gsvPoint.Y(); output << ")";
      }
      
      // print (xmin,ymin),(xmax,ymax) of VOC features
      for (unsigned int j = 0; j < images[i].vocFeatures.size(); j++) {
	output << " (";
	output << images[i].vocFeatures[j]->XMin(); output << ",";
	output << images[i].vocFeatures[j]->YMin(); output << "),(";
	output << images[i].vocFeatures[j]->XMax(); output << ",";
	output << images[i].vocFeatures[j]->YMax(); output << ")";
      }
      
      // print a new line after every image
      output << "\n";
    }
  }
  
  // close the output stream
  output.close();
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

int commandLineError() {
  // print how to use the program	
  fprintf(stderr, "Usage: correlateImages [-v] [-input_scene] <input_scene.gsv> [-input_gsv_correspondences] <gsv_correspondences.txt> [-input_voc_features] <voc_features.txt> -run int -segment int [-output_data <output_data.txt>]\n");
  fprintf(stderr, "\t-v: print statistics\n");
  fprintf(stderr, "\t-input_scene: .gsv scene file\n");
  fprintf(stderr, "\t-input_gsv_correspondences: .txt correspondence file from gsvalign\n");
  fprintf(stderr, "\t-input_voc_features: .txt voc feature file\n");
  fprintf(stderr, "\t-run: run number of correspondences\n");
  fprintf(stderr, "\t-segment: segment number of correspondences\n");
  fprintf(stderr, "\t-output_data: optional argument that outputs image data");
  
  // return FAILURE status
  return 0;
}

int parseArgs(int argc, char **argv) {
  // parse arguments
  argv++; argc--;
  while (argc > 0) {
    if (*argv[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-near_hit_range")) { argv++; argc--; nearHitRange = atof(*argv); }
      else if (!strcmp(*argv, "-input_scene")) { argv++; argc--; input_scene_name = *argv; }
      else if (!strcmp(*argv, "-input_gsv_correspondences")) { argv++; argc--; input_gsv_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-input_voc_features")) { 
	argv++; argc--;
	while (*argv[0] != '-') {
	  voc_files.push_back(*argv);
	  argv++; argc--;
	}
	argv--; argc++;
      }
      else if (!strcmp(*argv, "-output_data")) { argv++; argc--; output_data_name = *argv; }
      else if (!strcmp(*argv, "-run")) { argv++; argc--; run = atoi(*argv); }
      else if (!strcmp(*argv, "-segment")) { argv++; argc--; segment = atoi(*argv); }
      else if (!strcmp(*argv, "-usage")) { return commandLineError(); }
      else { fprintf(stderr, "Unrecognized command line argument: %s\n", *argv); return commandLineError(); }
    } 
    else {
      if (!input_scene_name) { input_scene_name = *argv; }
      else if (!input_gsv_correspondences_name) { input_gsv_correspondences_name = *argv; }
      else { voc_files.push_back(*argv); }
    }
    argv++; argc--;
  }
  
  // need an input
  if (!input_scene_name) { fprintf(stderr, "Need an input scene\n"); return commandLineError(); }
  // need a gsv correspondence file
  if (!input_gsv_correspondences_name) { fprintf(stderr, "Need a GSV correspondence file\n"); return commandLineError(); }
  // need a voc features file
  if (voc_files.size() == 0) { fprintf(stderr, "Need a VOC features file\n"); return commandLineError(); }
  // need to specify both run and segment numbers
  if (run == -1 || segment == -1) { fprintf(stderr, "Need to specify a run and segment number\n"); return commandLineError(); }
  
  // return OK status
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Main 
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // parse the command line arguments
  if (!parseArgs(argc, argv)) exit(-1);
  
  // read in the scene
  CorrelateImage *correlation = new CorrelateImage(run, segment, input_scene_name, input_gsv_correspondences_name);
  
  // get all precision and recall data for all voc files
  for (unsigned int i = 0; i < voc_files.size(); i++) {
    char *input_voc_features_name = voc_files[i];
    if (!correlation->AddVOCFeatures(input_voc_features_name)) { fprintf(stderr, "Failed to read %s\n", input_voc_features_name); exit(-1);}
    // get and print image recall
    double imageRecall = correlation->Recall();
    double imagePrecision = correlation->Precision();
    if (output_data_name) correlation->OutputData();
  }
  
  for (unsigned int i = 0; i < correlation->NImages(); i++) {
    if (correlation->Images(i).vocFeatures.size() != 0)
      printf("%d\n", i);
  }
	
  //findBBox(correlation, 847);
  double probability = MatchingBBoxes(correlation, 847, 838);
  findBBox(correlation, 847, 838);
  
  // return OK status
  return 1;
}
