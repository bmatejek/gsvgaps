// Include files 
#include "GSV/GSV.h"
#include "R2Shapes/R2Shapes.h"
#include <fglut/fglut.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

// Program arguments
static char *input_scene_name = NULL;
static GSVScene *scene = NULL;
static int print_verbose = 0;
static char *object = NULL;

// GLUT variables 
static int GLUTwindow = 0;
static int GLUTmodifiers = 0;
static int GLUTwindow_width = 484;
static int GLUTwindow_height = 648;

// Grid viewing variables
static R2Grid *image = NULL;
static R2Box image_window = R2null_box;
static R2Point image_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval image_range(0, 0);
static int ordered_index = 0;
static int GSV_IMAGE_WIDTH = 1936;
static int GSV_IMAGE_HEIGHT = 2592;

class OrderedPrediction {
public:
	OrderedPrediction(R2Box box, R2Box matchedBox, double score, bool hit, char filename[]);

public:
	R2Box box;
	R2Box matchedBox;
	double score;
	bool hit;
	char filename[4096];
};

OrderedPrediction::OrderedPrediction(R2Box box, R2Box matchedBox, double score, bool hit, char input_filename[]) :
box(box),
matchedBox(matchedBox),
score(score),
hit(hit)
{
	strcpy(filename, input_filename);
}

class TruePrediction {
public:
	TruePrediction(R2Box truth, int index);

public:
	R2Box truth;
	int index;
};

TruePrediction::TruePrediction(R2Box truth, int index) :
truth(truth),
index(index)
{
}

bool PredictionSort(OrderedPrediction *a, OrderedPrediction *b) {
	return a->score > b->score;
}

vector<OrderedPrediction *> predictions = vector<OrderedPrediction *>();


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *ReadScene(const char *filename) {
	// Start statistics
	RNTime start_time;
	start_time.Read();

	// Check if should read points (if .gsv, will be read as needed)
	RNBoolean read_points = (strstr(filename, ".gsv")) ? FALSE : TRUE;

	// Allocate google scene
	GSVScene *scene = new GSVScene();
	if (!scene) {
		fprintf(stderr, "Unable to allocate scene\n");
		return NULL;
	}

	// Read scene 
	if (!scene->ReadFile(filename, read_points)) {
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

static bool GridExists() {
	if (ifstream(predictions[ordered_index]->filename)) return true;
	else return false;
}

int ReadGrid(void) {
	if (!GridExists()) {
		fprintf(stderr, "Failed to read grid: %s\n", predictions[ordered_index]->filename);
		exit(-1);
	}

	// Allocate grid
	delete image;
	image = new R2Grid();
	image->Read(predictions[ordered_index]->filename);

	// return SUCCESS
	return 1;
}

void GLUTDrawText(const R2Point& p, const char *s) {
	// Draw text string s and position p
	glRasterPos2d(p[0], p[1]);
	while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}



void GLUTStop(void) {
	// Destroy window 
	glutDestroyWindow(GLUTwindow);

	// Exit
	exit(0);
}

void GLUTRedraw(void) {
	// Check grid
	if (!image) return;

	// Clear window 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set projection matrix
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(image_window.XMin(), image_window.XMax(), image_window.YMin(), image_window.YMax());

	// Set model view matrix
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// Draw grid values
	int xmin = (image_window.XMin() > 1) ? image_window.XMin() : 1;
	int ymin = (image_window.YMin() > 1) ? image_window.YMin() : 1;
	int xmax = (image_window.XMax() + 1 < image->XResolution() - 1) ? image_window.XMax() + 1 : image->XResolution() - 1;
	int ymax = (image_window.YMax() + 1 < image->YResolution() - 1) ? image_window.YMax() + 1 : image->YResolution() - 1;
	for (int j = ymin; j <= ymax; j++) {
		glBegin(GL_TRIANGLE_STRIP);
		for (int i = xmin; i <= xmax; i++) {
			for (int k = -1; k <= 0; k++) {
				RNScalar value = image->GridValue(i, j + k);
				RNLoadRgb(RNRgb(value, value, value));
				glVertex2i(i, j + k);
			}
		}
		glEnd();
	}

	// draw the lines for this bounding box
	R2Box box = predictions[ordered_index]->box;
	double xminp = box.XMin();
	double yminp = GSV_IMAGE_HEIGHT - box.YMax();
	double xmaxp = box.XMax();
	double ymaxp = GSV_IMAGE_HEIGHT - box.YMin();

	glBegin(GL_LINES);
	RNLoadRgb(RNRgb(1.0, 0.0, 0.0));
	glVertex3f(xminp, yminp, 0.0);
	glVertex3f(xmaxp, yminp, 0.0);

	glVertex3f(xmaxp, yminp, 0.0);
	glVertex3f(xmaxp, ymaxp, 0.0);

	glVertex3f(xmaxp, ymaxp, 0.0);
	glVertex3f(xminp, ymaxp, 0.0);

	glVertex3f(xminp, ymaxp, 0.0);
	glVertex3f(xminp, yminp, 0.0);
	glEnd();

	
	// draw the line for the matching box
	if (predictions[ordered_index]->hit) {
		R2Box matchingBox = predictions[ordered_index]->matchedBox;
		double xminp = matchingBox.XMin();
		double yminp = GSV_IMAGE_HEIGHT - matchingBox.YMax();
		double xmaxp = matchingBox.XMax();
		double ymaxp = GSV_IMAGE_HEIGHT - matchingBox.YMin();

		glBegin(GL_LINES);
		RNLoadRgb(RNRgb(0.0, 0.0, 1.0));
		glVertex3f(xminp, yminp, 0.0);
		glVertex3f(xmaxp, yminp, 0.0);

		glVertex3f(xmaxp, yminp, 0.0);
		glVertex3f(xmaxp, ymaxp, 0.0);

		glVertex3f(xmaxp, ymaxp, 0.0);
		glVertex3f(xminp, ymaxp, 0.0);

		glVertex3f(xminp, ymaxp, 0.0);
		glVertex3f(xminp, yminp, 0.0);
		glEnd();
	}

	// Reset projection matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	// Reset model view matrix
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	// Swap buffers 
	glutSwapBuffers();
}

void GLUTResize(int w, int h) {
	// Resize window
	glViewport(0, 0, w, h);

	// Remember window size 
	GLUTwindow_width = w;
	GLUTwindow_height = h;

	// Update selected grid window
	if (image) {
		RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
		RNScalar grid_aspect = (double)image->XResolution() / (double)image->YResolution();
		R2Point origin = image->GridBox().Centroid();
		R2Vector diagonal = image->GridBox().Max() - origin;
		diagonal[0] *= window_aspect / grid_aspect;
		image_window = R2Box(origin - diagonal, origin + diagonal);
	}

	// Redraw
	glutPostRedisplay();
}


void GLUTSpecial(int key, int x, int y) {
	// Process keyboard button event 
	switch (key) {
		case GLUT_KEY_UP: {
			ordered_index++;
			ordered_index = ordered_index % predictions.size();
			ReadGrid();
			break;
		}
		case GLUT_KEY_DOWN: {
			ordered_index--;
			if (ordered_index < 0) ordered_index = 0;
			ReadGrid();
			break;
		}
	}

	// Remember modifiers 
	GLUTmodifiers = glutGetModifiers();

	// Redraw
	glutPostRedisplay();
}

void GLUTKeyboard(unsigned char key, int x, int y) {
	// Process keyboard button event 
	switch (key) {			
		case 27: // ESCAPE
			GLUTStop();
			break;
	}

	// Remember modifiers 
	GLUTmodifiers = glutGetModifiers();

	// Redraw
	glutPostRedisplay();
}

void GLUTInit(int *argc, char **argv) {
	// Open window 
	glutInit(argc, argv);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
	GLUTwindow = glutCreateWindow("boxview");

	// Initialize background color 
	glClearColor(00, 0.0, 0.0, 1.0);

	// Initialize GLUT callback functions 
	glutDisplayFunc(GLUTRedraw);
	glutReshapeFunc(GLUTResize);
	glutKeyboardFunc(GLUTKeyboard);
	glutSpecialFunc(GLUTSpecial);
}

void GLUTMainLoop(void) {

	// Run main loop -- never returns 
	glutMainLoop();
}

static int Usage(void) {
	fprintf(stderr, "orderview input_scene_name [object]\n");
	return 0;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) input_scene_name = *argv;
			else if (!object) object = *argv;
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
			argv++; argc--;
		}
	}

	if (!input_scene_name) return Usage();
	if (!object) return Usage();

	// Return OK status 
	return 1;
}

static int OrderImages(void) {
	// go through every image
	for (int ir = 0; ir < scene->NRuns(); ir++) {
		GSVRun *run = scene->Run(ir);
		for (int is = 0, panorama_number = 0; is < run->NSegments(); is++) {
			GSVSegment *segment = run->Segment(is);
			for (int ip = 0; ip < segment->NPanoramas(); ++ip, ++panorama_number) {
				GSVPanorama *panorama = segment->Panorama(ip);
				for (int ii = 0; ii < panorama->NImages(); ii++) {
					
					// read in the projected truth for this image
					vector<TruePrediction *> proj_boxes = vector<TruePrediction *>();

					char projFilename[4096];
					sprintf(projFilename, "proj_truth/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

					ifstream truth;
					truth.open(projFilename);
					if (!truth.is_open()) { fprintf(stderr, "Failed to read %s\n", projFilename); return 0; }

					string line;
					while (getline(truth, line)) {
						vector<string> coordinates = splitString(line);
						if (coordinates.size() != 5) { fprintf(stderr, "Failed to read %s\n", projFilename); return 0; }
						double xmin = atof(coordinates[0].c_str());
						double ymin = atof(coordinates[1].c_str());
						double xmax = atof(coordinates[2].c_str());
						double ymax = atof(coordinates[3].c_str());

						int index = atoi(coordinates[4].c_str());

						R2Box box = R2Box(xmin, ymin, xmax, ymax);
						proj_boxes.push_back(new TruePrediction(box, index));
					}
					
					// get grid name
					char gridFilename[4096];
					sprintf(gridFilename, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run->Name(), is, panorama_number, ii);

					// read all of the predictions
					char  predFilename[4096];
					sprintf(predFilename, "voc_predictions/%s/%s/%02d_%06d_%02d_UndistortedImage.txt", object, run->Name(), is, panorama_number, ii);

					ifstream pred;
					pred.open(predFilename);
					if (!pred.is_open()) { fprintf(stderr, "Failed to open %s\n", predFilename); }

					while (getline(pred, line)) {
						vector<string> coordinates = splitString(line);
						if (coordinates.size() != 6) { fprintf(stderr, "Failed to read %s\n", predFilename); return 0; }

						double xmin = atof(coordinates[0].c_str());
						double ymin = atof(coordinates[1].c_str());
						double xmax = atof(coordinates[2].c_str());
						double ymax = atof(coordinates[3].c_str());

						R2Box box = R2Box(xmin, ymin, xmax, ymax);
						R2Box matchedBox = R2Box();
						double score = atof(coordinates[5].c_str());
						bool match = false;

						// go through every true bounding box in this image to see if there is a match
						for (unsigned int it = 0; it < proj_boxes.size(); it++) {
							R2Box union_truth = proj_boxes[it]->truth;
							R2Box intersection_truth = proj_boxes[it]->truth;

							union_truth.Union(box);
							intersection_truth.Intersect(box);

							// see if the intersection over union is greater than 0.5
							if (intersection_truth.Area() / union_truth.Area() > 0.10) {
								match = true;
								matchedBox = proj_boxes[it]->truth;
								break;
							}
						}

						OrderedPrediction *ordered = new OrderedPrediction(box, matchedBox, score, match, gridFilename);
						predictions.push_back(ordered);
					}
				}
			}
		}
	}

	sort(predictions.begin(), predictions.end(), PredictionSort);

	// return success
	return 1;
}

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// read the GSV Scene
	scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);

	// get ordered list of objects
	if (!OrderImages()) exit(-1);

	// Read grid
	if (!ReadGrid()) exit(-1);

	// Initialize GLUT
	GLUTInit(&argc, argv);

	// Run GLUT interface
	GLUTMainLoop();

	// Return success 
	return 0;
}