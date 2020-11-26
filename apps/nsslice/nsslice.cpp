// Source file for the mesh viewer program



// Include files 
#include "R3Graphics/R3Graphics.h"
#include "R3GridArray.h"
#include "R2Shapes/R2Shapes.h"
#include <fglut/fglut.h>



// Program variables

static char *input_name;
static int print_verbose = 0;
static int print_debug = 0;
static int animate = 0;

// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// Grid viewing variables

static R3GridArray *grids = NULL;
static R3Grid *grid = NULL;
static R2Grid *selected_slice = NULL;
static int selected_slice_index = -1;
static R2Box selected_slice_window = R2null_box;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval selected_slice_range(0, 0);
static int color_type = 1; // 0=gray, 1=red-green-blue



static void
SelectGrid(int index)
{
    // Check index
    if (index < 0) index = 0;
    if (index > grid->ZResolution() - 1) index = grid->ZResolution() - 1;
    if (index == selected_slice_index) return;
    
    R2Grid *slice = grid->Slice(2, index);
    // Set window title
    char title[4096];
    sprintf(title, "Slice %d", index);
    glutSetWindowTitle(title);

    // Update display variables
    if (!selected_slice || (selected_slice->XResolution() != slice->XResolution()) || (selected_slice->YResolution() != slice->YResolution()) ||  (0 && (selected_slice->WorldToGridTransformation() != slice->WorldToGridTransformation()))) {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)grid->XResolution() / (double)grid->YResolution();
        R2Point origin = slice->GridBox().Centroid();
        R2Vector diagonal = slice->GridBox().Max() - origin;
        diagonal[0] *= window_aspect / grid_aspect;
        selected_slice_window = R2Box(origin - diagonal, origin + diagonal);
        selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
    }

    // Update min and max values
    selected_slice_range = grid->Range();

    // Update selected grid 
    selected_slice_index = index;
    selected_slice = slice;
}



static RNRgb
Color(RNScalar value)
{
    // Check for unknown value
    if (value == R2_GRID_UNKNOWN_VALUE) {
        if (color_type == 0) return RNRgb(1, 0.5, 0);
        else return RNblack_rgb;
    }

    // Normalize color
    RNScalar value_min = selected_slice_range.Min();
    RNScalar value_width = selected_slice_range.Diameter();
    RNScalar value_scale = (value_width > 0) ? 1.0 / value_width : 1.0;
    RNScalar normalized_value = value_scale * (value - value_min);

    // Compute color
    RNRgb c(0, 0, 0);
    if (color_type == 0) {
        c[0] = normalized_value;
        c[1] = normalized_value;
        c[2] = normalized_value;
    }
    else {
        if (normalized_value < 0.5) {
            c[0] = 1 - 2 * normalized_value;
            c[1] = 2 * normalized_value;
        }
        else {
            c[1] = 1 - 2 * (normalized_value - 0.5);
            c[2] = 2 * (normalized_value - 0.5);
        }
    }

    // Return color
    return c;
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
    if (!selected_slice) return;

    // Clear window 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set projection matrix
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(selected_slice_window.XMin(), selected_slice_window.XMax(), selected_slice_window.YMin(), selected_slice_window.YMax());

    // Set model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Draw grid values
    int xmin = (selected_slice_window.XMin() > 1) ? selected_slice_window.XMin() : 1;
    int ymin = (selected_slice_window.YMin() > 1) ? selected_slice_window.YMin() : 1;
    int xmax = (selected_slice_window.XMax() + 1 < selected_slice->XResolution() - 1) ? selected_slice_window.XMax() + 1 : selected_slice->XResolution() - 1;
    int ymax = (selected_slice_window.YMax() + 1 < selected_slice->YResolution() - 1) ? selected_slice_window.YMax() + 1 : selected_slice->YResolution() - 1;
    for (int j = ymin; j <= ymax; j++) {
        glBegin(GL_TRIANGLE_STRIP);
        for (int i = xmin; i <= xmax; i++) {
            for (int k = -1; k <= 0; k++) {
                RNScalar value = selected_slice->GridValue(i, j + k);
                RNRgb color = Color(value);
                RNLoadRgb(color);
                glVertex2i(i, j + k);
            }
        }
        glEnd();
    }

    // Draw value at selected grid position
    if ((selected_slice_position.X() != RN_UNKNOWN) && (selected_slice_position.Y() != RN_UNKNOWN)) {
        int ix = (int)(selected_slice_position.X() + 0.5);
        int iy = (int)(selected_slice_position.Y() + 0.5);
        RNScalar value = selected_slice->GridValue(ix, iy);
        char buffer[1024];
        if (value != R2_GRID_UNKNOWN_VALUE) sprintf(buffer, "%d %d : %g", ix, iy, value);
        else sprintf(buffer, "%d %d : %s", ix, iy, "Unknown");
        RNRgb color = Color(value);
        RNRgb complement = RNwhite_rgb - color;
        RNLoadRgb(RNmagenta_rgb);
        R2Box(selected_slice_position - 0.5 * R2ones_vector, selected_slice_position + 0.5 * R2ones_vector);
        GLUTDrawText(selected_slice_position + 2 * R2ones_vector, buffer);
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

void timer(int v) {
    if (animate) {
        if (selected_slice_index < grid->ZResolution() - 2) {
            SelectGrid(selected_slice_index + 1);
            GLUTRedraw();
        }
        else {
            Sleep(2000);
            SelectGrid(0);
            GLUTRedraw();
        }
    }
    glutTimerFunc(50, timer, v);
}


void GLUTResize(int w, int h)
{
    // Resize window
    glViewport(0, 0, w, h);

    // Remember window size 
    GLUTwindow_width = w;
    GLUTwindow_height = h;

    // Update selected grid window
    if (selected_slice) {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)selected_slice->XResolution() / (double)selected_slice->YResolution();
        R2Point origin = selected_slice->GridBox().Centroid();
        R2Vector diagonal = selected_slice->GridBox().Max() - origin;
        diagonal[0] *= window_aspect / grid_aspect;
        selected_slice_window = R2Box(origin - diagonal, origin + diagonal);
    }

    // Redraw
    glutPostRedisplay();
}



void GLUTMotion(int x, int y)
{
    // Invert y coordinate
    y = GLUTwindow_height - y;

    // Compute mouse movement
    int dx = x - GLUTmouse[0];
    int dy = y - GLUTmouse[1];

    // View manipulation
    if (selected_slice) {
        if (GLUTbutton[0]) {
            // Query
            RNScalar px = x * selected_slice_window.XLength() / (double)GLUTwindow_width + selected_slice_window.XMin();
            RNScalar py = y * selected_slice_window.YLength() / (double)GLUTwindow_height + selected_slice_window.YMin();
            selected_slice_position.Reset(px, py);
            glutPostRedisplay();
        }
        else if (GLUTbutton[1]) {
            // Zoom
            RNScalar scale_factor = 1;
            scale_factor *= 1.0 - (double)dx / (double)GLUTwindow_width;
            scale_factor *= 1.0 - (double)dy / (double)GLUTwindow_height;
            scale_factor *= scale_factor;
            selected_slice_window.Inflate(scale_factor);
            glutPostRedisplay();
        }
        else if (GLUTbutton[2]) {
            // Pan
            RNScalar tx = -dx * selected_slice_window.XLength() / (double)GLUTwindow_width;
            RNScalar ty = -dy * selected_slice_window.YLength() / (double)GLUTwindow_height;
            selected_slice_window.Translate(R2Vector(tx, ty));
            glutPostRedisplay();
        }
    }

    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
    // Invert y coordinate
    y = GLUTwindow_height - y;

    // Process mouse button event
    if (button == 0) {
        if (state == GLUT_DOWN) {
            // Query
            RNScalar px = x * selected_slice_window.XLength() / (double)GLUTwindow_width + selected_slice_window.XMin();
            RNScalar py = y * selected_slice_window.YLength() / (double)GLUTwindow_height + selected_slice_window.YMin();
            selected_slice_position.Reset(px, py);
            glutPostRedisplay();
        }
        else {
            selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
            glutPostRedisplay();
        }
    }
    else if ((button == 3) || (button == 4)) {
        if (state == GLUT_DOWN) {
            // Zoom with wheel
            RNScalar scale_factor = (button == 3) ? 0.9 : 1.1;
            selected_slice_window.Inflate(scale_factor);
            glutPostRedisplay();
        }
    }

    // Remember button state 
    int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
    GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

    // Remember modifiers 
    GLUTmodifiers = glutGetModifiers();

    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // Redraw
    glutPostRedisplay();
}



void GLUTSpecial(int key, int x, int y)
{
    // Invert y coordinate
    y = GLUTwindow_height - y;

    // Process keyboard button event 
    switch (key) {
    case GLUT_KEY_PAGE_UP:
    case GLUT_KEY_PAGE_DOWN:
        if (selected_slice) {
            int shift = 0;
            if (key == GLUT_KEY_PAGE_UP) shift = -1;
            else if (key == GLUT_KEY_PAGE_DOWN) shift = 1;
            SelectGrid(selected_slice_index + shift);
            glutPostRedisplay();
        }
        break;

    case GLUT_KEY_LEFT:
    case GLUT_KEY_RIGHT:
    case GLUT_KEY_UP:
    case GLUT_KEY_DOWN:
        if (selected_slice) {
            RNScalar selected_slice_minimum = selected_slice->Minimum();
            RNScalar selected_slice_maximum = selected_slice->Maximum();
            RNScalar selected_slice_radius = 0.5 * (selected_slice_maximum - selected_slice_minimum);
            RNScalar center = selected_slice_range.Mid();
            RNScalar radius = selected_slice_range.Radius();
            if (key == GLUT_KEY_LEFT) radius *= 0.95;
            else if (key == GLUT_KEY_RIGHT) radius *= 1.05;
            else if (key == GLUT_KEY_UP) center += 0.01 * selected_slice_radius;
            else if (key == GLUT_KEY_DOWN) center -= 0.01 * selected_slice_radius;
            if (radius > selected_slice_radius) radius = selected_slice_radius;
            if (center - radius < selected_slice_minimum) center = selected_slice_minimum + radius;
            if (center + radius > selected_slice_maximum) center = selected_slice_maximum - radius;
            RNScalar min_value = center - radius;
            RNScalar max_value = center + radius;
            selected_slice_range.Reset(min_value, max_value);
            glutPostRedisplay();
        }
        break;
    }

    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // Remember modifiers 
    GLUTmodifiers = glutGetModifiers();

    // Redraw
    glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
    // Invert y coordinate
    y = GLUTwindow_height - y;

    // Process keyboard button event 
    switch (key) {
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
        SelectGrid(key - '1');
        break;

    case 'C':
    case 'c':
        color_type = ((color_type + 1) % 2);
        break;

    case 'A':
    case 'a':
        animate = 1 - animate;
        break;

    case 27: // ESCAPE
        GLUTStop();
        break;
    }

    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // Remember modifiers 
    GLUTmodifiers = glutGetModifiers();

    // Redraw
    glutPostRedisplay();
}




#if 0

void GLUTIdle(void)
{
    // Set current window
    if (glutGetWindow() != GLUTwindow)
        glutSetWindow(GLUTwindow);

    // Redraw
    glutPostRedisplay();
}

#endif



void GLUTInit(int *argc, char **argv)
{
    // Set window dimensions
    // TODO make this cleaner
    GLUTwindow_width = 800;
    GLUTwindow_height = 800;
    if (!grids->NGrids()) {
        /// TODO Make cleaner
        grid = grids->Grid(0);
        RNScalar aspect = (RNScalar)grid->YResolution() / (RNScalar)grid->XResolution();
        GLUTwindow_height = aspect * GLUTwindow_width;
    }
    // Open window 
    glutInit(argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
    GLUTwindow = glutCreateWindow("nsslice");

    // Initialize background color 
    glClearColor(00, 0.0, 0.0, 1.0);
    // Initialize GLUT callback functions 
    glutDisplayFunc(GLUTRedraw);
    glutReshapeFunc(GLUTResize);
    glutKeyboardFunc(GLUTKeyboard);
    glutSpecialFunc(GLUTSpecial);
    glutMouseFunc(GLUTMouse);
    glutMotionFunc(GLUTMotion);
    glutTimerFunc(100, timer, 0);
}



void GLUTMainLoop(void)
{
    // Select first grid
    SelectGrid(0);
    // Run main loop -- never returns 
    glutMainLoop();
}


static R3GridArray *
ReadGridArray(char *input_name)
{
    // Start statistics
    RNTime start_time;
    start_time.Read();

    // Allocate a grid array
    R3GridArray *grids = new R3GridArray();
    if (!grids) {
        fprintf(stderr, "Unable to allocate grid\n");
        return NULL;
    }

    // Read grids
    int status = grids->ReadFile(input_name);
    if (!status) {
        fprintf(stderr, "Unable to read grid file %s\n", input_name);
        return NULL;
    }

    // Check grids
    if (grids->NGrids() == 0) {
        fprintf(stderr, "Grid file is empty %s\n", input_name);
        delete grids;
        return NULL;
    }

    // Print statistics
    if (print_verbose) {
        printf("Read grid ...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  # Grids = %d\n", grids->NGrids());
        fflush(stdout);
    }

    // Print info for each grid
    if (print_debug) {
        for (int i = 0; i < grids->NGrids(); i++) {
            R3Grid *grid = grids->Grid(i);
            RNInterval grid_range = grid->Range();
            printf("  Grid %d:\n", i);
            printf("    Resolution = %d %d %d\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
            printf("    Spacing = %g\n", grid->GridToWorldScaleFactor());
            printf("    Cardinality = %d\n", grid->Cardinality());
            printf("    Volume = %g\n", grid->Volume());
            printf("    Minimum = %g\n", grid_range.Min());
            printf("    Maximum = %g\n", grid_range.Max());
            printf("    L1Norm = %g\n", grid->L1Norm());
            printf("    L2Norm = %g\n", grid->L2Norm());
            fflush(stdout);
        }
    }

    // Return grids
    return grids;
}


////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
    // Parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) print_debug = 1;
            else if (!strcmp(*argv, "-animate")) animate = 1;
            else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        }
        else {
            if (!input_name) input_name = *argv;
            else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        }
        argv++; argc--;
    }

    // Check filenames
    if (!input_name) {
        fprintf(stderr, "Usage: nsview gridfile [text] [options]\n");
        return 0;
    }

    // Return OK status 
    return 1;
}


int main(int argc, char **argv)
{
    // Parse program arguments
    if (!ParseArgs(argc, argv)) exit(-1);

    // Read grid
    grids = ReadGridArray(input_name);
    if (!grids) exit(-1);

    // read in the first grid
    grid = grids->Grid(0);

    // Initialize GLUT
    GLUTInit(&argc, argv);

    // Run GLUT interface
    GLUTMainLoop();

    // Return success 
    return 0;
}








