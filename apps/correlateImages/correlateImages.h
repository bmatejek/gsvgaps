#ifndef CORRELATE_IMAGES_H
#define CORRELATE_IMAGES_H

////////////////////////////////////////////////////////////////////////
// Include Files
////////////////////////////////////////////////////////////////////////

#include "vocFeature.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

struct GSVLidarFeature {
	GSVFeature *lidar;
	bool matched;
	bool nearlyMatched;
};

struct LidarImagePair {
	GSVLidarFeature *gsvLidarFeature;
	GSVFeature *gsvImageFeature;
};

struct Image {
	std::vector<VOCFeature *> vocFeatures;
	std::vector<GSVFeature *> gsvFeatures;
	std::vector<LidarImagePair *> lidarImagePairs;
	GSVImage *image;
};

class CorrelateImage {
public:
	// constructor
	CorrelateImage(int, int, char *, char *);
	
	// Accessors/Manipulators
	GSVPoseOptimization *Optimization(void);
	Image Images(int i);
	const int NImages(void) const;
	double Recall(void);
	double Precision(void);
	int AddVOCFeatures(char *);
	void OutputData(void);
	
private:
	// instance variables
	int run;
	int segment;
	GSVPoseOptimization *optimization;
	VOCFeatureList *vocFeatures;
	std::vector<Image> images;
	std::vector<struct GSVLidarFeature *> gsvLidarFeatures;
	int numberOfImages;
};

#endif // CORRELATE_IMAGES_H
