// Source file for the mesh viewer program



// Include files 

#include "R2Shapes/R2Shapes.h"
#include <fglut/fglut.h>
#include <string>
#include <vector>
#include <fstream>

using namespace std;



// Program variables

static char grid_name[4096];
static int print_verbose = 0;



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 648;
static int GLUTwindow_width = 484;
static int GLUTmodifiers = 0;



// Grid viewing variables

static R2Grid *image = NULL;
static R2Box image_window = R2null_box;
static R2Point image_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval image_range(0, 0);

// bounding box object
class BBox {
public:
	BBox(vector<R2Point> points, string object, int index, bool difficult, bool truncated);

public:
	vector<R2Point> points;
	string object;
	int index;
	bool difficult;
	bool truncated;
};

BBox::BBox(vector<R2Point> points, string object, int index, bool difficult, bool truncated) :
points(points),
object(object),
index(index),
difficult(difficult),
truncated(truncated)
{
}

static char run[128];
static int segment_index = 0;
static int panorama_index = 0;
static int image_index = 0;
static vector<R2Point> current = vector<R2Point>();
static vector<BBox> boxes = vector<BBox>();
static string objects[5] = { "traffic_light", "stop_sign", "street_light", "one_way_sign", "signal_box" };
static RNRgb colors[5] = { RNRgb(0.0, 1.0, 1.0), RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), RNRgb(1.0, 1.0, 0.0) };
static int object_index = 0;
static bool difficult = false;
static bool truncated = false;

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



int ReadGrid(void) {
	current.clear();
	boxes.clear();

	// Allocate grid
	delete image;
	image = new R2Grid();
	assert(image);

	// Read grid from file
	if (!image->Read(grid_name)) {
		delete image;

		// return failure
		return 0;
	}

	// return SUCCESS
	return 1;
}

void GLUTDrawText(const R2Point& p, const char *s)
{
	// Draw text string s and position p
	glRasterPos2d(p[0], p[1]);
	while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}



void GLUTStop(void)
{
	// Destroy window 
	glutDestroyWindow(GLUTwindow);

	// Exit
	exit(0);
}



void GLUTRedraw(void)
{
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
				glColor3f(value, value, value);
				glVertex2i(i, j + k);
			}
		}
		glEnd();
	}

	// draw in all of the bboxes
	
	// draw the current bbox
	if (current.size() > 2) {
		glBegin(GL_LINES);
		RNLoadRgb(colors[object_index]);
		glVertex2f(4 * current[0].X(), 2592 - 4 * current[0].Y());
		for (unsigned int c = 1; c < current.size(); c++) {
			glVertex2f(4 * current[c].X(), 2592 - 4 * current[c].Y());

			if (c != current.size() - 1)
				glVertex2f(4 * current[c].X(), 2592 - 4 * current[c].Y());
		}
		glEnd();
	}
	else if (current.size() == 2) {
		glBegin(GL_LINES);
		RNLoadRgb(colors[object_index]);
		glVertex2f(4 * current[0].X(), 2592 - 4 * current[0].Y());
		glVertex2f(4 * current[1].X(), 2592 - 4 * current[1].Y());
		glEnd();
	}

	// go through all of the boxes
	for (unsigned int b = 0; b < boxes.size(); b++) {
		glBegin(GL_LINES);
		RNLoadRgb(colors[boxes[b].index]);
		vector<R2Point> points = boxes[b].points;
		glVertex2f(4 * points[0].X(), 2592 - 4 * points[0].Y());
		for (unsigned int p = 0; p < points.size(); p++) {
			glVertex2f(4 * points[p].X(), 2592 - 4 * points[p].Y());

			glVertex2f(4 * points[p].X(), 2592 - 4 * points[p].Y());
		}
		glVertex2f(4 * points[0].X(), 2592 - 4 * points[0].Y());
		glEnd();
	}

	//// Reset projection matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	//// Reset model view matrix
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	
	// Swap buffers 
	glutSwapBuffers();
}



void GLUTResize(int w, int h)
{
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

static void SaveBoundingBoxes(void) {
	char save_name[4096];
	sprintf(save_name, "truth_boxes/%s/%02d_%06d_%02d_UndistortedImage.txt", run, segment_index, panorama_index, image_index);

	if (boxes.size() == 0) return;
	printf("Saving bounding boxes in %s\n", save_name);

	// actually save the bounding boxes
	ofstream ofs;
	ofs.open(save_name, std::ofstream::trunc);
	if (!ofs.is_open()) {
		fprintf(stderr, "Unable to open %s\n", save_name);
		return;
	}

	for (unsigned int b = 0; b < boxes.size(); b++) {
		vector<R2Point> points = boxes[b].points;
		ofs << boxes[b].object;
		ofs << '\n';
		ofs << boxes[b].difficult;
		ofs << '\n';
		ofs << boxes[b].truncated;
		ofs << '\n';
		for (unsigned int p = 0; p < points.size(); p++) {
			char output[10];
			sprintf(output, "%d,%d\n", (int)round(4 * points[p].X()), (int)round(4 * points[p].Y()));
			ofs << output;
		}
	}

	ofs.close();
}

void GLUTSpecial(int key, int x, int y)
{

	// Process keyboard button event 
	switch (key) {
	case GLUT_KEY_LEFT:
		SaveBoundingBoxes();
		image_index--;
		if (image_index < 0) image_index = 7;
		sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		ReadGrid();
		difficult = false;
		truncated = false;
		printf("%s\n", grid_name);
		break;
	case GLUT_KEY_RIGHT:
		SaveBoundingBoxes();
		image_index++;
		if (image_index == 8) image_index = 0;
		sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		ReadGrid();
		difficult = false;
		truncated = false;
		printf("%s\n", grid_name);
		break;
	case GLUT_KEY_UP:
		SaveBoundingBoxes();
		panorama_index++;
		sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		while (!ReadGrid()) {
			panorama_index--;
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		}
		difficult = false;
		truncated = false;
		printf("%s\n", grid_name);
		break;
	case GLUT_KEY_DOWN:
		SaveBoundingBoxes();
		panorama_index--;
		sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		while (!ReadGrid()) {
			panorama_index++;
			sprintf(grid_name, "gsv_data/%s/%02d_%06d_%02d_UndistortedImage.jpg", run, segment_index, panorama_index, image_index);
		}
		difficult = false;
		truncated = false;
		printf("%s\n", grid_name);
		break;
	}

	// Remember modifiers 
	GLUTmodifiers = glutGetModifiers();

	// Redraw
	glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{

	// Process keyboard button event 
	switch (key) {			
		case '1':
			object_index = 0;
			printf("%s\n", objects[object_index].c_str());
			break;
		case '2':
			object_index = 1;
			printf("%s\n", objects[object_index].c_str());
			break;
		case '3':
			object_index = 2;
			printf("%s\n", objects[object_index].c_str());
			break;
		case '4':
			object_index = 3;
			printf("%s\n", objects[object_index].c_str());
			break;
		case '5':
			object_index = 4;
			printf("%s\n", objects[object_index].c_str());
			break;

		case 'D':
		case 'd':
			difficult = !difficult;
			printf("Difficult? %d\n", difficult);
			break;

		case 'T':
		case 't':
			truncated = !truncated;
			printf("Truncated? %d\n", truncated);
			break;

		// save all of the bounding boxes
		case 19:
			SaveBoundingBoxes();
			break;

		case 26: // Ctrl + Z
			if (current.size() == 0 && boxes.size() != 0) boxes.pop_back();
			if (current.size() != 0) current.pop_back();
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

void GLUTMouse(int button, int state, int x, int y) {
	if (state == GLUT_UP) {
		switch (button) {
			case GLUT_LEFT_BUTTON: {
				current.push_back(R2Point(x, y));
				break;
			}
			case GLUT_RIGHT_BUTTON: {
				if (current.size() == 0) break;
				boxes.push_back(BBox(current, objects[object_index], object_index, difficult, truncated));
				current.clear();
				difficult = false;
				truncated = false;
				break;
			}
		}

		// Redraw
		glutPostRedisplay();
	}
}

void GLUTInit(int *argc, char **argv)
{
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
	glutMouseFunc(GLUTMouse);
}



void GLUTMainLoop(void)
{

	// Run main loop -- never returns 
	glutMainLoop();
}


static int
ParseArgs(int argc, char **argv)
{
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
			argv++; argc--;
		}
		else {
			strcpy(grid_name, *argv);
			argv++; argc--;
		}
	}

	// Return OK status 
	return 1;
}

static int ValidateImage(void) {
	vector<string> fileparts = splitString(grid_name, '/');
	vector<string> indices = splitString(fileparts[2].c_str(), '_');

	strcpy(run, fileparts[1].c_str());
	segment_index = atoi(indices[0].c_str());
	panorama_index = atoi(indices[1].c_str());
	image_index = atoi(indices[2].c_str());

	// return SUCCESS
	return 1;
}

int main(int argc, char **argv)
{
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	// validate the image
	if (!ValidateImage()) exit(-1);

	// Read grid
	if (!ReadGrid()) exit(-1);

	// Initialize GLUT
	GLUTInit(&argc, argv);

	// Run GLUT interface
	GLUTMainLoop();

	// Return success 
	return 0;
}








