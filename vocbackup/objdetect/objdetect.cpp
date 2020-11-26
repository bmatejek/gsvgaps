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

using namespace std;

////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// Input/output
static char *input_scene_name = NULL;
static int print_verbose = 0;
static int scanline = 0;
static double resolution = 1.0;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

static char *object = NULL;
static char *output_filename = NULL;
static GSVScene *scene = NULL;
static R3Grid *output = NULL;
static R2Grid **image_probabilities = NULL;
static double max_distance = 20.0;

// Useful constants
static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;
static double low_threshold = -0.5;
static double high_threshold = 1.0;


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

static int CreateOutput(void) {
	// determine the size of the gsv input
	R3Box bbox = scene->Run(0)->BBox();
	for (int ir = 1; ir < scene->NRuns(); ir++) {
		bbox.Union(scene->Run(ir)->BBox());
	}

	int x = ceil(resolution * bbox.XLength());
	int y = ceil(resolution * bbox.YLength());
	int z = ceil(resolution * bbox.ZLength());

	output = new R3Grid(x, y, z, bbox);
	if (!output) return 0;

	return 1;
}

static int SaveProbabilityGrid(void) {
	RNTime start_time;
	start_time.Read();

	output->WriteGridFile(output_filename);

	if (print_verbose) {
		printf("Saved object predictions to %s ...\n", output_filename);
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
	}

	// return OK status
	return 1;
}

// split the line by delim and return the vector
static vector<string> splitString(string line, char delim = ',') {
	vector<string> arguments;
	int length = line.length();
	int string_start = 0;

	// split the string by delim
	for (int i = 0; i < length; i++) {
		if (line.at(i) == delim) {
			int string_end = i;
			if (string_start != string_end) {
				arguments.push_back(line.substr(string_start, string_end - string_start));
			}
			string_start = i + 1;
		}
	}

	// add the last substring
	if (string_start != length)
		arguments.push_back(line.substr(string_start));

	// return the number of arguments found
	return arguments;
}

static double ConvertScore(double score) {
	if (score < low_threshold) return 0.0;
	else if (score > high_threshold) return 1.0;
	else return (score - low_threshold) / (high_threshold - low_threshold);
}

static int CreateImageProbabilities(void) {
	RNTime start_time;
	start_time.Read();

	image_probabilities = new R2Grid *[scene->NImages()];

	// read in all of the voc probabilities (FOR NOW)
	int i = 0;
	for (int ir = 0; ir < scene->NRuns(); ++ir) {
		GSVRun *run = scene->Run(ir);
		int panorama_number = 0;
		for (int is = 0; is < run->NSegments(); ++is) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ++ii, ++i) {
					image_probabilities[i] = new R2Grid(GSV_IMAGE_WIDTH, GSV_IMAGE_HEIGHT);
					if (!image_probabilities[i]) return 0;

					char probabilityFilename[128];
					sprintf(probabilityFilename, "voc_predictions/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

					ifstream probabilities;
					probabilities.open(probabilityFilename);
					if (!probabilities.is_open()) { fprintf(stderr, "Failed to open %s\n", probabilityFilename); return 0; }

					string line;
					while (getline(probabilities, line)) {
						vector<string> prediction = splitString(line);
						
						// get coordinates of box
						double xmin = atof(prediction[0].c_str());
						double ymax = GSV_IMAGE_HEIGHT - atof(prediction[1].c_str());
						double xmax = atof(prediction[2].c_str());
						double ymin = GSV_IMAGE_HEIGHT - atof(prediction[3].c_str());
						double score = atof(prediction[5].c_str());

						double probability = ConvertScore(score);

						// add in this probability 
						for (int x = floor(xmin); x <= ceil(xmax); x++) {
							for (int y = floor(ymin); y <= ceil(ymax); y++) {
								double old_probability = image_probabilities[i]->GridValue(x, y);
								if (old_probability < probability) {
									image_probabilities[i]->SetGridValue(x, y, probability);
								}
							}
						}

					}
					probabilities.close();
				}
			}
		}
	}

	if (print_verbose) {
		printf("Read all probability distributions ... \n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
	}

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Object Detection Functions
////////////////////////////////////////////////////////////////////////

static int DetectObjects(void) {
	int i = 0; int n = 0;
	GSVScan *scan = NULL;

	// find the scanline we want
	for (int ir = 0; ir < scene->NRuns(); ++ir) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0; is < run->NSegments(); ++is) {
			GSVSegment *segment = run->Segment(is);
			for (int ia = 0; ia < segment->NScans(); ++ia, ++i) {
				if (i == scanline) {
					scan = segment->Scan(ia);
					if (!scan->ReadPoints()) {
						fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), scan->RunIndex(), scan->SegmentIndex());
						return 0;
					}

					for (int il = 0; il < scan->NScanlines(); ++il) {
						GSVScanline *scanline = scan->Scanline(il);
						n += scanline->NPoints();
					}

				}
			}
		}
	}

	if (scan == NULL) {
		fprintf(stderr, "Unable to find this scanline: %d\n", scanline);
		return 0;
	}

	// get all of the points for this scanline
	for (int il = 0, i = 0; il < scan->NScanlines(); ++il) {
		GSVScanline *scanline = scan->Scanline(il);
		for (int ip = 0; ip < scanline->NPoints(); ++ip, ++i) {
			// print out the percent of program that  is done
			printf("%0.6f%%\n", 100 * (double)i / n);

			R3Point world = scanline->PointPosition(il);

			// probability of this added grid resolution
			double add_probability = 0.0;

			// for every pair of images
			int count = 0;
			for (int ir = 0; ir < scene->NRuns(); ++ir) {
				GSVRun *run = scene->Run(ir);
				int panorama_number = 0;
				for (int is = 0; is < run->NSegments(); ++is) {
					GSVSegment *segment = run->Segment(is);
					for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
						GSVPanorama *panorama = segment->Panorama(ip);
						for (int ii = 0; ii < panorama->NImages(); ++ii, ++count) {
							GSVImage *image = panorama->Image(ii);

							// only consider if this image is within max_distance of object
							R3Point viewpoint = image->Pose().Viewpoint();
							if (R3Distance(viewpoint, world) > max_distance) continue;

							R2Point point = image->UndistortedPosition(world);
							double probabiliity = image_probabilities[count]->GridValue(point.X(), point.Y());

							// get information for second image
							int countp = 0;
							for (int irp = 0; irp < scene->NRuns(); ++irp) {
								GSVRun *runp = scene->Run(irp);
								int panorama_numberp = 0;
								for (int isp = 0; isp < run->NSegments(); ++isp) {
									GSVSegment *segmentp = runp->Segment(isp);
									for (int ipp = 0; ipp < segmentp->NPanoramas(); ++ipp, ++panorama_numberp) {
										GSVPanorama *panoramap = segmentp->Panorama(ipp);
										for (int iip = 0; iip < panoramap->NImages(); ++iip, ++countp) {
											// only want to consider half of the images
											if (countp >= count) goto nested_loop;

											GSVImage *imagep = panoramap->Image(iip);

											R2Point pointp = imagep->UndistortedPosition(world);

											// only consider if this image is within max_distancce of object
											R3Point viewpointp = imagep->Pose().Viewpoint();
											if (R3Distance(viewpointp, world) > max_distance) continue;

											double probabilityp = image_probabilities[countp]->GridValue(pointp.X(), pointp.Y());
											add_probability += probabiliity * probabilityp;

											// make sure that the two images are never the same
											if (ir == irp && is == isp && ip == ipp && ii == iip) {
												fprintf(stderr, "Image double counted %d %d %d %d\n", ir, is, ip, ii);
											}
										}
									}
								}
							}
						nested_loop:;
						}
					}
				}
			}
		}
	}

	// release all the scan  points
	scan->ReleasePoints();

	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage(void) {
	fprintf(stderr, "Usage: objdetect input_scenefile output_filename object [-v] [-resolution <resolution>] [-max_distance <max_distance>]\n");

	// return failure status
	return 0;
}

static int ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
	  else if (!strcmp(*argv, "-resolution")) { argv++; argc--; resolution = atof(*argv); }
	  else if (!strcmp(*argv, "-object")) { argv++; argc--; object = *argv; }
	  else if (!strcmp(*argv, "-max_distance")) { argv++; argc--; max_distance = atof(*argv); }
	  else if (!strcmp(*argv, "-scanline")) { argv++; argc--; scanline = atoi(*argv); }
	  else if (!strcmp(*argv, "-usage")) { return Usage(); }
      else { 
        fprintf(stderr, "Invalid program argument: %s\n", *argv); 
		return Usage();
      }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
	  else if (!output_filename) output_filename = *argv;
	  else if (!object) object = *argv;
	  else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) return Usage();

  // make sure there is an object
  if (!object) return Usage();

  // make sure there is an output filename
  if (!output_filename) return Usage();

  // Return OK status 
  return 1;
}

////////////////////////////////////////////////////////////////////////
// Evaluation Functions
////////////////////////////////////////////////////////////////////////

struct score_t {
	double score;
	bool hit;
};

bool sort_function(score_t a, score_t b) {
	return a.score < b.score;
}

static int EvaluatePredictions(void) {
	vector<R3Point> true_locations = vector<R3Point>();

	// get all truth data points
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);

		char truthFilename[128];
		sprintf(truthFilename, "aligned_truth/%s/%s/truth.txt", object, run->Name());
		
		ifstream truth;
		truth.open(truthFilename);

		if (!truth.is_open()) { fprintf(stderr, "Failed to read truth filename: %s\n", truthFilename); return 0; }

		string line;
		while (getline(truth, line)) {
			vector<string> location = splitString(line);

			double x = atof(location[0].c_str());
			double y = atof(location[1].c_str());
			double z = atof(location[2].c_str());

			true_locations.push_back(R3Point(x, y, z));
		}

		truth.close();
	}

	// used to exclude nearzero  scores
	double epsilon = 10e-6;
	vector<score_t> scores = vector<score_t>();
	for (int x = 0; x < output->XResolution(); ++x) {
		for (int y = 0; y < output->YResolution(); ++y) {
			for (int z = 0; z < output->ZResolution(); ++z) {
				double score = output->GridValue(x, y, z);
				
				R3Point world = output->WorldPosition(x, y, z);

				bool hit = false;
				// see if this score matches any of the true locations within a distance of the resolution?
				for (unsigned int i = 0; i < true_locations.size(); i++) {
					if (R3Distance(world, true_locations[i]) < 3) {
						hit = true;
					}
				}

				if (score > epsilon) {
					score_t st;	st.score = score; st.hit = hit;
					scores.push_back(st);
					printf("%f %f %f\n", world.X(), world.Y(), world.Z());
				}
			}
		}
	}

	sort(scores.begin(), scores.end(), sort_function);
	//printf("%d\n", scores.size());
	for (int i = 0; i < scores.size(); i++) {
		//printf("%0.4f %d\n", scores[i].score, scores[i].hit);
	}

	// return OK status
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	RNTime start_time;
	start_time.Read();

	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// Read scene
	scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);

	// read in all of the grids
	if (!CreateImageProbabilities()) exit(-1);
	
	// create the probability output
	if (!CreateOutput()) exit(-1);

	// detect the objects
	if (!DetectObjects()) exit(-1);

	// save the output grid
	if (!SaveProbabilityGrid()) exit(-1);

	// evaluate how well the predictions did
	if (!EvaluatePredictions()) exit(-1);

	if (print_verbose) {
		printf("Created probability distribution ... \n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
	}

	// Return success 
	return 0;
}


