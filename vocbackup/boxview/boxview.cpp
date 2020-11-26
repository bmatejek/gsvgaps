// Include files 
#include "R2Shapes/R2Shapes.h"
#include <fglut/fglut.h>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

// Program arguments
static int print_verbose = 0;
static char *object = NULL;

// Program variables
static char grid_name[4096];

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
static double low_threshold = -0.5;

static char run[128];
static int segment_index = 0;
static int panorama_index = 0;
static int image_index = 0;

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
	if (ifstream(grid_name)) return true;
	else return false;
}

int ReadGrid(void) {
	if (!GridExists()) {
		fprintf(stderr, "Cannot advance image to %s\n", grid_name);
		return 0;
	}
	
	// Allocate grid
	delete image;
	image = new R2Grid();
	image->Read(grid_name);

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

	// read in prediction boxes
	vector<string> fileparts = splitString(grid_name, '/');
	vector<string> filename = splitString(fileparts[2], '.');
	string prediction_file = "voc_predictions/" + string(object) + "/" + fileparts[1] + "/" + filename[0] + ".txt";
	ifstream predictions;
	predictions.open(prediction_file.c_str());
	if (!predictions.is_open()) {
		fprintf(stderr, "Failed to open %s\n", prediction_file.c_str());
	}

	// read every bbox in the file
	string line;
	while (getline(predictions, line)) {
		vector<string> prediction = splitString(line);

		// get coordinates of box
		int xmin = floor(atof(prediction[0].c_str()));
		int ymax = ceil(2592 - atof(prediction[1].c_str()));
		int xmax = ceil(atof(prediction[2].c_str()));
		int ymin = floor(2592 - atof(prediction[3].c_str()));
		double score = atof(prediction[5].c_str());

		if (score > low_threshold) {
			// draw the lines for this bounding box
			glBegin(GL_LINES);
			RNLoadRgb(RNRgb(1.0, 0.0, 0.0));
			glVertex3f(xmin, ymin, 0.0);
			glVertex3f(xmax, ymin, 0.0);

			glVertex3f(xmax, ymin, 0.0);
			glVertex3f(xmax, ymax, 0.0);

			glVertex3f(xmax, ymax, 0.0);
			glVertex3f(xmin, ymax, 0.0);

			glVertex3f(xmin, ymax, 0.0);
			glVertex3f(xmin, ymin, 0.0);
			glEnd();
		}
	}
	
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
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			ReadGrid();
			break;
		case GLUT_KEY_RIGHT:
			image_index++;
			if (image_index == 8) image_index = 0;
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			ReadGrid();
			break;
		case GLUT_KEY_UP:
			panorama_index++;
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			while (!ReadGrid()) {
				panorama_index--;
				sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			}
			break;
		case GLUT_KEY_DOWN:
			panorama_index--;
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			while (!ReadGrid()) {
				panorama_index++;
				sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
			}
			break;
	}

	// Remember modifiers 
	GLUTmodifiers = glutGetModifiers();

	// Redraw
	glutPostRedisplay();
}

void GLUTKeyboard(unsigned char key, int x, int y) {
	// Process keyboard button event 
	switch (key) {			
		case 'S':
		case 's':
			low_threshold -= 0.1;
			break;

		case 'W':
		case 'w':
			low_threshold += 0.1;
			break;

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
	fprintf(stderr, "boxview [object] [first_image]\n");
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
			if (!object) object = *argv;
			strcpy(grid_name, *argv);
			argv++; argc--;
		}
	}

	if (!object) return Usage();

	// Return OK status 
	return 1;
}

static int ValidateImage(void) {
	vector<string> fileparts = splitString(grid_name, '/');
	if (fileparts.size() != 3) return 0;

	vector<string> indices = splitString(fileparts[2].c_str(), '_');
	if (indices.size() != 4) return 0;

	strcpy(run, fileparts[1].c_str());
	segment_index = atoi(indices[0].c_str());
	panorama_index = atoi(indices[1].c_str());
	image_index = atoi(indices[2].c_str());
	// return SUCCESS
	return 1;
}

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// validate the image
	if (!ValidateImage()) { fprintf(stderr, "Unable to read %s\n", grid_name); exit(-1); }

	// Read grid
	if (!ReadGrid()) exit(-1);

	// Initialize GLUT
	GLUTInit(&argc, argv);

	// Run GLUT interface
	GLUTMainLoop();

	// Return success 
	return 0;
}