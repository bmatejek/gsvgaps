////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "helper.h"

using namespace std;

#define SUCCESS 1
#define FAILURE 0
#define GSV_WIDTH 1936
#define GSV_HEIGHT 2592

////////////////////////////////////////////////////////////////////////
// Command line arguments
////////////////////////////////////////////////////////////////////////

// command line arguments
int print_verbose  = 0;
int image_maxima = 1;
char *input_scene_name = NULL;
GSVScene *scene = NULL;
double max_distance = 40;

// converts the VOC prediction score into a probability
static double ConvertProbability(double value) {
	if (value < 0.0)
		return 0.0;
	else if (value > 1.0)
		return 1.0;
	else
		return value;
}

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

// write grid to filename
static void WriteGrid(const char *filename, R3Grid *grid) {
	RNTime start_time;
	start_time.Read();

	// write the grid
	if (!grid->WriteGridFile(filename)) {
		fprintf(stderr, "Unable to write %s\n", filename);
		return;
	}

	// print statistics
	if (print_verbose) {
		printf("Write grid to %s ...\n", filename);
		printf("\tTime = %.2f seconds\n", start_time.Elapsed());
	}
}

// read in the GSVScene from filename
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

////////////////////////////////////////////////////////////////////////
// Generate the Grids
////////////////////////////////////////////////////////////////////////

static void RasterizeWorldSpan(R3Grid *grid, R3Ray ray, double low_estimate, double high_estimate, double value) {
	R3Point start = grid->GridPosition(ray.Point(low_estimate));
	R3Point end = grid->GridPosition(ray.Point(high_estimate));

	const double *pre1 = start.Coords();
	const double *pre2 = end.Coords();

	int p1[3] = { (int) round(pre1[0]), (int) round(pre1[1]), (int) round(pre1[2]) };
	int p2[3] = { (int) round(pre2[0]), (int) round(pre2[1]), (int) round(pre2[2]) };

	// Get some convenient variables
	int d[3],p[3],dd[3],s[3];
	for (int i = 0; i < 3; i++) {
		d[i]= p2[i] - p1[i];
		if(d[i]<0){
			dd[i] = -d[i];
			s[i] = -1;
		}
		else{
			dd[i] = d[i];
			s[i] = 1;
		}
		p[i] = p1[i];
	}

	// Choose dimensions
	int i1=0;
	if(dd[1]>dd[i1]){i1=1;}
	if(dd[2]>dd[i1]){i1=2;}
	int i2=(i1+1)%3;
	int i3=(i1+2)%3;

	// Check span extent
	if(dd[i1]==0){
		// Span is a point - rasterize it
		if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
			double distance = R3Distance(ray.Start(), R3Point(p[0], p[1], p[2]));
			grid->AddGridValue(p[0], p[1], p[2], value * distance * distance);
		}
	}
	else {
		// Step along span
		int off[3] = { 0, 0, 0 };
		for (int i = 0; i <= dd[i1]; i++) {
			if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
				double distance  = R3Distance(ray.Start(), R3Point(p[0], p[1], p[2]));
				grid->AddGridValue(p[0], p[1], p[2], value * distance * distance);
			}
			off[i2]+=dd[i2];
			off[i3]+=dd[i3];
			p[i1]+=s[i1];
			p[i2]+=s[i2]*off[i2]/dd[i1];
			p[i3]+=s[i3]*off[i3]/dd[i1];
			off[i2]%=dd[i1];
			off[i3]%=dd[i1];
		}
	}
}

static void RasterizeWorldSpan(R3Grid *grid, R3Ray ray, double max_distance, double value) {
	R3Point *intersect_point = new R3Point();
	R3Vector normal = R3Vector();
	RNScalar dist = RNScalar();

	double max_diagonal = 1.2 * sqrt((grid->WorldBox().XMax() - grid->WorldBox().XMin()) * (grid->WorldBox().XMax() - grid->WorldBox().XMin()) + (grid->WorldBox().YMax() - grid->WorldBox().YMin()) * (grid->WorldBox().YMax() - grid->WorldBox().YMin()) + (grid->WorldBox().ZMax() - grid->WorldBox().ZMin()) * (grid->WorldBox().ZMax() - grid->WorldBox().ZMin()));

	R3Ray secondary = R3Ray(ray.Point(-1 * max_diagonal), ray.Vector());
	R3Intersects(grid->WorldBox(), secondary, intersect_point, &normal, &dist);

	R3Point start = grid->GridPosition(ray.Start());
	R3Point end;
	if (R3Distance(grid->GridPosition(*intersect_point), start) > max_distance)
		end = grid->GridPosition(ray.Point(max_distance));
	else
		end = grid->GridPosition(*intersect_point);

	const double *pre1 = start.Coords();
	const double *pre2 = end.Coords();

	int p1[3] = { (int) round(pre1[0]), (int) round(pre1[1]), (int) round(pre1[2]) };
	int p2[3] = { (int) round(pre2[0]), (int) round(pre2[1]), (int) round(pre2[2]) };

	// Get some convenient variables
	int d[3],p[3],dd[3],s[3];
	for (int i = 0; i < 3; i++) {
		d[i]= p2[i] - p1[i];
		if(d[i]<0){
			dd[i] = -d[i];
			s[i] = -1;
		}
		else{
			dd[i] = d[i];
			s[i] = 1;
		}
		p[i] = p1[i];
	}

	// Choose dimensions
	int i1=0;
	if(dd[1]>dd[i1]){i1=1;}
	if(dd[2]>dd[i1]){i1=2;}
	int i2=(i1+1)%3;
	int i3=(i1+2)%3;

	// Check span extent
	if(dd[i1]==0){
		// Span is a point - rasterize it
		if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
			double distance = R3Distance(start, R3Point(p[0], p[1], p[2]));
			grid->AddGridValue(p[0], p[1], p[2], value * distance * distance);
		}
	}
	else {
		// Step along span
		int off[3] = { 0, 0, 0 };

		for (int i = 0; i <= dd[i1]; i++) {
			if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
				double distance  = R3Distance(start, R3Point(p[0], p[1], p[2]));
				grid->AddGridValue(p[0], p[1], p[2], value * distance * distance);
			}
			off[i2]+=dd[i2];
			off[i3]+=dd[i3];
			p[i1]+=s[i1];
			p[i2]+=s[i2]*off[i2]/dd[i1];
			p[i3]+=s[i3]*off[i3]/dd[i1];
			off[i2]%=dd[i1];
			off[i3]%=dd[i1];
		}
	}
	delete intersect_point;
}

static void SecondaryVoting(R3Grid *hough, R3Grid *grid, R3Ray ray, double max_distance, double value) {
	R3Point *intersect_point = new R3Point();
	R3Vector normal = R3Vector();
	RNScalar dist = RNScalar();

	double max_diagonal = 1.2 * sqrt((grid->WorldBox().XMax() - grid->WorldBox().XMin()) * (grid->WorldBox().XMax() - grid->WorldBox().XMin()) + (grid->WorldBox().YMax() - grid->WorldBox().YMin()) * (grid->WorldBox().YMax() - grid->WorldBox().YMin()) + (grid->WorldBox().ZMax() - grid->WorldBox().ZMin()) * (grid->WorldBox().ZMax() - grid->WorldBox().ZMin()));

	R3Ray secondary = R3Ray(ray.Point(-1 * max_diagonal), ray.Vector());
	R3Intersects(grid->WorldBox(), secondary, intersect_point, &normal, &dist);

	R3Point start = grid->GridPosition(ray.Start());
	R3Point end;
	if (R3Distance(grid->GridPosition(*intersect_point), start) > max_distance)
		end = grid->GridPosition(ray.Point(max_distance));
	else
		end = grid->GridPosition(*intersect_point);

	const double *pre1 = start.Coords();
	const double *pre2 = end.Coords();

	int p1[3] = { (int) round(pre1[0]), (int) round(pre1[1]), (int) round(pre1[2]) };
	int p2[3] = { (int) round(pre2[0]), (int) round(pre2[1]), (int) round(pre2[2]) };

	// Get some convenient variables
	int d[3],p[3],dd[3],s[3];
	for (int i = 0; i < 3; i++) {
		d[i]= p2[i] - p1[i];
		if(d[i]<0){
			dd[i] = -d[i];
			s[i] = -1;
		}
		else{
			dd[i] = d[i];
			s[i] = 1;
		}
		p[i] = p1[i];
	}

	// Choose dimensions
	int i1=0;
	if(dd[1]>dd[i1]){i1=1;}
	if(dd[2]>dd[i1]){i1=2;}
	int i2=(i1+1)%3;
	int i3=(i1+2)%3;

	double max_value = 0.0;
	R3Point max_point = R3Point();

	// Check span extent
	if(dd[i1]==0){
		// Span is a point - rasterize it
		if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
			if (max_value < grid->GridValue(p[0], p[1], p[2])) {
				max_value = grid->GridValue(p[0], p[1], p[2]);
				max_point = R3Point(p[0], p[1], p[2]);
			}
		}
	}
	else {
		// Step along span
		int off[3] = { 0, 0, 0 };

		for (int i = 0; i <= dd[i1]; i++) {
			if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
				if (max_value < grid->GridValue(p[0], p[1], p[2])) {
					max_value = grid->GridValue(p[0], p[1], p[2]);
					max_point = R3Point(p[0], p[1], p[2]);
				}
			}
			off[i2]+=dd[i2];
			off[i3]+=dd[i3];
			p[i1]+=s[i1];
			p[i2]+=s[i2]*off[i2]/dd[i1];
			p[i3]+=s[i3]*off[i3]/dd[i1];
			off[i2]%=dd[i1];
			off[i3]%=dd[i1];
		}
	}
	double distance = R3Distance(max_point, ray.Start());
	double add_value = value * distance * distance;
	hough->AddGridValue(max_point.X(), max_point.Y(), max_point.Z(), add_value);
	delete intersect_point;
}

static void SecondaryVoting(R3Grid *hough, R3Grid *grid, R3Ray ray, double low_estimate, double high_estimate, double value) {
	R3Point start = grid->GridPosition(ray.Point(low_estimate));
	R3Point end = grid->GridPosition(ray.Point(high_estimate));

	const double *pre1 = start.Coords();
	const double *pre2 = end.Coords();

	int p1[3] = { (int) round(pre1[0]), (int) round(pre1[1]), (int) round(pre1[2]) };
	int p2[3] = { (int) round(pre2[0]), (int) round(pre2[1]), (int) round(pre2[2]) };

	// Get some convenient variables
	int d[3],p[3],dd[3],s[3];
	for (int i = 0; i < 3; i++) {
		d[i]= p2[i] - p1[i];
		if(d[i]<0){
			dd[i] = -d[i];
			s[i] = -1;
		}
		else{
			dd[i] = d[i];
			s[i] = 1;
		}
		p[i] = p1[i];
	}

	// Choose dimensions
	int i1=0;
	if(dd[1]>dd[i1]){i1=1;}
	if(dd[2]>dd[i1]){i1=2;}
	int i2=(i1+1)%3;
	int i3=(i1+2)%3;

	double max_value = 0.0;
	R3Point max_point = R3Point();

	// Check span extent
	if(dd[i1]==0){
		// Span is a point - rasterize it
		if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
			if (max_value < grid->GridValue(p[0], p[1], p[2])) {
				max_value = grid->GridValue(p[0], p[1], p[2]);
				max_point = R3Point(p[0], p[1], p[2]);
			}
		}
	}
	else {
		// Step along span
		int off[3] = { 0, 0, 0 };

		for (int i = 0; i <= dd[i1]; i++) {
			if (((p[0] >= 0) && (p[0] < grid->Resolution(0))) && ((p[1] >= 0) && (p[1] < grid->Resolution(1))) && ((p[2] >= 0) && (p[2] < grid->Resolution(2)))) {
				if (max_value < grid->GridValue(p[0], p[1], p[2])) {
					max_value = grid->GridValue(p[0], p[1], p[2]);
					max_point = R3Point(p[0], p[1], p[2]);
				}
			}
			off[i2]+=dd[i2];
			off[i3]+=dd[i3];
			p[i1]+=s[i1];
			p[i2]+=s[i2]*off[i2]/dd[i1];
			p[i3]+=s[i3]*off[i3]/dd[i1];
			off[i2]%=dd[i1];
			off[i3]%=dd[i1];
		}
	}
	double distance = R3Distance(max_point, ray.Start());
	double add_value = value * distance * distance;
	hough->AddGridValue(max_point.X(), max_point.Y(), max_point.Z(), add_value);
}

static int ForwardHoughVote(void) {
	int scale_factor = 16;
	bool range = false;

	// start statistics
	RNTime start_time;
	start_time.Read();

	// go through all runs in this scene
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);

		// get existing bounding box
		R3Box bbox = run->BBox();

		double scaling_factor = 4.0;
		double max_distance = scaling_factor * 30.0;
		int xresolution = (int) (scaling_factor * round(bbox.XMax() - bbox.XMin()));
		int yresolution = (int) (scaling_factor * round(bbox.YMax() - bbox.YMin()));
		int zresolution = (int) (scaling_factor * round(bbox.ZMax() - bbox.ZMin()));

		// print the grid resolution
		if (print_verbose)
			printf("Grid Resolution: %d %d %d\n", xresolution, yresolution, zresolution);

		// create new R3Grid for preliminary and secondary voting
		R3Grid *grid = new R3Grid(xresolution, yresolution, zresolution, bbox);
		R3Grid *hough = new R3Grid(xresolution, yresolution, zresolution, bbox);

		// go through each image in this run
		int panorama_number = 0;
		for (int is = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
					GSVImage *image = panorama->Image(ii);

					char probabilityFilename[128];
					snprintf(probabilityFilename, 128, "voc_probabilities/traffic_light/%s/%02d_%06d_%02d_UndistortedImage.grd", run->Name(), is, panorama_number, ii);


					R2Grid *values = new R2Grid();
					values->Read(probabilityFilename);

					printf("Preliminary Voting: %02d_%06d_%02d\n", is, panorama_number, ii);
					for (int x = 0; x < GSV_WIDTH / scale_factor; x++) {
						for (int y = 0; y < GSV_HEIGHT / scale_factor; y++) {
							R3Ray ray = image->RayThroughUndistortedPosition(R2Point(x * scale_factor, y * scale_factor));

							double value = ConvertProbability(values->GridValue(x, y));

							if (value != 0.0) {
								if (range) {
									/* TODO ADD IN WAY TO READ THE LOW AND THE HIGH ESTIMATE */
									double low_estimate = 0;
									double high_estimate = max_distance;
									RasterizeWorldSpan(grid, ray, low_estimate, high_estimate, value);
								} else {
									RasterizeWorldSpan(grid, ray, max_distance, value);
								}
							}
						}
					}

					delete values;
				}
			}
		}

		// iterate through all of the images again to mask non maxima along the rays
		panorama_number = 0;
		for (int is = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages() - 1; ii++) {
					GSVImage *image = panorama->Image(ii);

					char probabilityFilename[128];
					snprintf(probabilityFilename, 128, "voc_probabilities/traffic_light/%s/%02d_%06d_%02d_UndistortedImage.grd", run->Name(), is, panorama_number, ii);

					R2Grid *values = new R2Grid();
					values->Read(probabilityFilename);

					printf("Masking Non Maxima: %02d_%06d_%02d\n", is, panorama_number, ii);
					for (int x = 0; x < GSV_WIDTH / scale_factor; x++) {
						for (int y = 0; y < GSV_HEIGHT / scale_factor; y++) {
							R3Ray ray = image->RayThroughUndistortedPosition(R2Point(x, y));

							double value = ConvertProbability(values->GridValue(x, y));

							if (value != 0.0) {
								if (range) {
									/* TODO ADD IN WAY TO READ THE LOW AND THE HIGH ESTIMATES */
									double low_estimate = 0;
									double high_estimate = max_distance;
									SecondaryVoting(hough, grid, ray, low_estimate, high_estimate, value);
								} else {
									SecondaryVoting(hough, grid, ray, max_distance, value);
								}
							}
						}
					}
				}
			}
		}
		WriteGrid("hough.grd", hough);
		delete hough;
		delete grid;
	}
	return SUCCESS;
}

static int BackwardHoughVote(void) {
	// images are downsampled by a factor of 16
	int scale_factor = 16;
	bool range = false;

	// start statistics
	RNTime start_time;
	start_time.Read();

	// go through all runs in this scene
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);

		// get existing bounding box
		R3Box bbox = run->BBox();

		// grid samples are 1 meter apart
		double scaling_factor = 1.0;
		double max_distance = scaling_factor * 30.0;
		int xresolution = (int) (scaling_factor * round(bbox.XMax() - bbox.XMin()));
		int yresolution = (int) (scaling_factor * round(bbox.YMax() - bbox.YMin()));
		int zresolution = (int) (scaling_factor * round(bbox.ZMax() - bbox.ZMin()));

		// print the grid resolution
		if (print_verbose)
			printf("Grid Resolution: %d %d %d\n", xresolution, yresolution, zresolution);

		// create new R3Grid for preliminary and secondary voting
		R3Grid *grid = new R3Grid(xresolution, yresolution, zresolution, bbox);
		int panorama_number = 0;
		for (int is = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ii++) {
					GSVImage *image = panorama->Image(ii);

					char probabilityFilename[128];
					snprintf(probabilityFilename, 128, "voc_probabilities/traffic_light/%s/%02d_%06d_%02d_UndistortedImage.grd", run->Name(), is, panorama_number, ii);

					R2Grid *values = new R2Grid();
					values->Read(probabilityFilename);

					// find all of the grid samples that are viewed by this image
					printf("%02d_%06d_%02d\n", is, panorama_number, ii);
					for (int x = 0; x < xresolution; x++) {
						for (int y = 0; y < yresolution; y++) {
							for (int z = 0; z < zresolution; z++) {
								R3Point world_position = grid->WorldPosition(x, y, z);
								if (R3Distance(image->Pose().Viewpoint(), world_position) < max_distance) {
									R2Point point = image->UndistortedPosition(world_position);
									if (point.X() >= 0 && point.X() < GSV_WIDTH && point.Y() >= 0 && point.Y() < GSV_HEIGHT) {
										// add in this grid value
										double value = ConvertProbability(values->GridValue(point.X() / scale_factor, point.Y() / scale_factor));
										grid->AddGridValue(x, y, z, value);
									}
								}			
							}
						}
					}
					delete values;
				}
			}
		}
		WriteGrid("backwards_hough.grd", grid);
		delete grid;
	}
	return SUCCESS;
}



////////////////////////////////////////////////////////////////////////
// Argument parsing functions
////////////////////////////////////////////////////////////////////////

// return how this program should be used
static int Usage(void) {
	fprintf(stderr, "hough.exe [-v] <scene.gsv>\n");
	return FAILURE;
}

// parse the command line arguments
static int ParseArgs(int argc, char **argv) {
// parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-scene")) { argv++; input_scene_name = *argv; }
			else { fprintf(stderr, "Invalid program argument: %s.\n", *argv); return FAILURE; }
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
		}
		argc--; argv++;
	}

	if (!input_scene_name) { fprintf(stderr, "Input scene required\n"); return FAILURE; }

	return SUCCESS;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	// start statistics
	RNTime start_time;
	start_time.Read();

	// parse function arguments
	if (!ParseArgs(argc, argv)) { return Usage(); }

	// read scene
	scene = ReadScene(input_scene_name);
	if (!scene) { fprintf(stderr, "Unable to read scene: %s\n", input_scene_name); return FAILURE; }

	if (false) {
		for (int ir = 0; ir < scene->NRuns(); ir++) {
			GSVRun *run = scene->Run(ir);
			int panorama_number = 0;
			for (int is = 0; is < run->NSegments(); is++) {
				GSVSegment *segment = run->Segment(is);
				for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
					GSVPanorama *panorama = segment->Panorama(ip);
					for (int ii = 0; ii < panorama->NImages(); ii++) {
						// create new R2Grid
						R2Grid *grid = new R2Grid(GSV_WIDTH / 16, GSV_HEIGHT / 16);
						for (int x = 0; x < GSV_WIDTH / 16; x++) {
							for (int y = 0; y < GSV_HEIGHT / 16; y++) {
								grid->SetGridValue(x, y, -2.0);
							}
						}

						char filename[128];
						snprintf(filename, 128, "voc_predictions/traffic_light/%s/%02d_%06d_%02d_UndistortedImage.txt", run->Name(), is, panorama_number, ii);
						printf("%s\n", filename);
						ifstream file(filename);
						string line;
						while (getline(file, line)) {
							vector<string> box = splitString(line, ',');
							int xmin = (int) round(atof(box[0].c_str())) / 16;
							int ymax = (int) (GSV_HEIGHT - round(atof(box[1].c_str()))) / 16;
							int xmax = (int) round(atof(box[2].c_str())) / 16;
							int ymin = (int) (GSV_HEIGHT - round(atof(box[3].c_str()))) / 16; 
							double value = atof(box[5].c_str());
							for (int x = xmin; x < xmax; x++) {
								for (int y = ymin; y < ymax; y++) {
									if (grid->GridValue(x,y) < value) {
										grid->SetGridValue(x, y, value);
									}
								}
							}
						}
						char printFile[128];
						snprintf(printFile, 128, "voc_probabilities/traffic_light/%s/%02d_%06d_%02d_UndistortedImage.grd", run->Name(), is, panorama_number, ii);
						grid->WriteGridFile(printFile);
						delete grid;
					}
				}
			}
		}
	}

	// call the hough voting algorithm
	//if (!ForwardHoughVote()) { fprintf(stderr, "Failed to cast hough votes\n"); return FAILURE; }

	// call the backwards hough voting algorithm
	if (!BackwardHoughVote()) { fprintf(stderr, "Failed to cast hough votes\n"); return FAILURE; }

	// print statistics
	if (print_verbose) {
		printf("Program Completed\n");
		printf("\tTime = %.2f seconds\n", start_time.Elapsed());
	}

	return SUCCESS;
}
