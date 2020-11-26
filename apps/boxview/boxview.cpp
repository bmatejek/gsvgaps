// Include files 
#include "R2Shapes/R2Shapes.h"
#include <fglut/fglut.h>
#include <GSV/GSV.h>
#include <string>
#include <vector>
#include <fstream>
#include <R3Shapes/R3Shapes.h>

using namespace std;

// Program arguments
static int print_verbose = 0;
static char *input_scene_name = NULL;
static GSVScene *scene = NULL;
static GSVRun *run = NULL;

// Global variables
int image_index = 0;
int segment_index = 0;
int panorama_index = 0;
int segment_panorama_index = 0;

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
static double threshold = -1.0;
static int GSV_IMAGE_HEIGHT = 2592;

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

static bool GridExists(char *gridFilename) {
	if (ifstream(gridFilename)) return true;
	else return false;
}

int ReadGrid(void) {
	char gridFilename[4096];
	sprintf(gridFilename, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run->Name(), segment_index, panorama_index, image_index);

	if (!GridExists(gridFilename)) {
		fprintf(stderr, "Failed to read grid: %s\n", gridFilename);
		exit(-1);
	}

	// Allocate grid
	delete image;
	image = new R2Grid();
	image->Read(gridFilename);

	// return SUCCESS
	return 1;
}


static void DrawTextRCNN(const R3Point& p, const char *s, void *font = GLUT_BITMAP_HELVETICA_10) {
	// Draw text string s and position p
	glRasterPos3d(p[0], p[1], p[2]);
	while (*s) glutBitmapCharacter(font, *(s++));
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

	// read in prediction boxes
	char filename[4096];
	sprintf(filename, "rcnn/%s/%02d_%06d_%02d_UndistortedImage.txt", run->Name(), segment_index, panorama_index, image_index);
	
	ifstream predictions;
	predictions.open(filename);
	if (!predictions.is_open()) { fprintf(stderr, "Failed to open %s\n", filename); }

	// read every bbox in the file
	string line;
	while (getline(predictions, line)) {
		vector<string> prediction = splitString(line);

		// get coordinates of box
		string object = prediction[0];
		double score = atof(prediction[1].c_str());
		if (score < threshold) continue;
		double xminp = atof(prediction[2].c_str());
		double yminp = GSV_IMAGE_HEIGHT - atof(prediction[5].c_str());
		double xmaxp = atof(prediction[4].c_str());
		double ymaxp = GSV_IMAGE_HEIGHT - atof(prediction[3].c_str());

		// draw the lines for this bounding box
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

		RNLoadRgb(RNRgb(0.0, 1.0, 0.0));
		DrawTextRCNN(R3Point(xminp, ymaxp, 0.0), object.c_str());
		char cScore[10];
		sprintf(cScore, "%0.2lf", score);
		DrawTextRCNN(R3Point(xminp, yminp, 0.0), cScore);
	}

	char cThreshold[128];
	sprintf(cThreshold, "Threshold: %0.2lf", threshold);
	DrawTextRCNN(R3Point(0.0, 0.95 * GSV_IMAGE_HEIGHT, 0.0), cThreshold, GLUT_BITMAP_HELVETICA_18);

	// close the file
	predictions.close();

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
		case GLUT_KEY_LEFT:
			image_index--;
			if (image_index < 0) image_index = 7;
			break;
		case GLUT_KEY_RIGHT:
			image_index++;
			if (image_index == 8) image_index = 0;
			break;
		case GLUT_KEY_UP:
			panorama_index++;
			segment_panorama_index++;
			if (run->Segment(segment_index)->NPanoramas() <= segment_panorama_index) {
				segment_index++;
				segment_panorama_index = 0;
				if (segment_index >= run->NSegments()) {
					segment_index--;
					panorama_index--;
					segment_panorama_index = run->Segment(segment_index)->NPanoramas();
				}
			}
			break;
		case GLUT_KEY_DOWN:
			panorama_index--;
			segment_panorama_index--;
			if (segment_panorama_index < 0) {
				segment_index--;
				if (segment_index < 0) {
					segment_index = 0;
					segment_panorama_index = 0;
					panorama_index = 0;
				}
				else {
					segment_panorama_index = run->Segment(segment_index)->NPanoramas();
				}
			}
			break;
		case GLUT_KEY_PAGE_UP:
			threshold += 0.1;
			break;
		case GLUT_KEY_PAGE_DOWN:
			threshold -= 0.1;
			if (threshold < -1.0) { threshold = -1.0; }
			break;
	}

	ReadGrid();

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
	fprintf(stderr, "boxview input_scene_name [-v]\n");
	return 0;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
			if (!input_scene_name) { input_scene_name = *argv; }
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// require input scene
	if (!input_scene_name) return Usage();

	// Return OK status 
	return 1;
}

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// read the GSV scene
	scene = ReadScene(input_scene_name);
	if (!scene) exit(-1);

	// only allow scenes with one run
	if (scene->NRuns() > 1) { fprintf(stderr, "Only implemented for scenes with 1 run\n"); exit(-1); }
	run = scene->Run(0);

	// read the first grid
	ReadGrid();

	// Initialize GLUT
	GLUTInit(&argc, argv);

	// Run GLUT interface
	GLUTMainLoop();

	// Return success 
	return 0;
}