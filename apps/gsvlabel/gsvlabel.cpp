// Source file for the GSV labeling program



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "object.h"
#include <fglut/fglut.h>



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static const char *input_scene_name = NULL;
static const char *input_labels_name = NULL;
static const char *input_segmentation_name = "Segmentation";
static const char *input_ground_truth_name = NULL;
static const char *input_keystroke_name = NULL;
static const char *output_keystroke_name = NULL;
static const char *output_scan_rectangle_name = NULL;
static const char *output_image_rectangle_name = NULL;
static RNArray<const char *> channel_names;
static const char *parameterization = "DA";
static int interactive = 0;
static int print_verbose = 0;



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

enum {
  NO_COLOR_SCHEME,
  RGB_COLOR_SCHEME,
  CHANNEL_COLOR_SCHEME,
  SIMILARITY_COLOR_SCHEME,
  SEGMENTATION_COLOR_SCHEME,
  LABEL_COLOR_SCHEME,
  PICK_COLOR_SCHEME,
  NUM_COLOR_SCHEMES
};

struct StreetSide {
public:
  StreetSide(GSVScan *scan);
  ~StreetSide(void);
  int Read(void);
  int Write(void);
  void Release(void);
  void Draw(int color_scheme = RGB_COLOR_SCHEME) const;
  void DrawImage(int color_scheme = RGB_COLOR_SCHEME) const;
  void DrawPoints(int color_scheme = RGB_COLOR_SCHEME) const;
  void DrawGrid(int color_scheme = RGB_COLOR_SCHEME) const;
  R3Viewer PointViewer(void) const; 
public:
  GSVScan *scan;
  ObjectParse *parse;
  R2Grid *scanline_grid;
  ObjectPoint **point_backpointers;
  ObjectSegmentation **segmentation_backpointers;
  RNArray<struct Channel *> channels;
  int index;
};

struct Channel {
public:
  Channel(StreetSide *streetside, const char *parameterization, const char *name, const char *extension = "grd");
  ~Channel(void);
  int Read(void);
  void Release(void);
public:
  StreetSide *streetside;  
  char *parameterization;
  char *name;
  char *extension;
  RNArray<R2Grid *> grids;
  RNScalar mean, stddev;
  RNScalar minimum, maximum;
  int width, height;
  int index;
};

struct Classifier {
public:
  Classifier(void);
  ~Classifier(void);
public:
  int Read(const char *filename);
  int Write(const char *filename) const;
  int Train(const RNArray<StreetSide *>& streetsides, const RNArray<const char *>& classifier_feature_names);
  int Predict(const RNArray<StreetSide *>& streetsides, const RNArray<const char *>& classifier_feature_names, const RNArray<const char *>& crf_feature_names);
};

struct Keystroke {
  Keystroke(StreetSide *streetside, const char *parameterization, const R2Point& p, int key) 
    : streetside(streetside), parameterization((parameterization) ? strdup(parameterization) : NULL), p(p), key(key) {};
  void Draw(void) const;
  StreetSide *streetside;
  char *parameterization;
  R2Point p;
  int key;
};

const int max_segmentation_colors = 12;
const RNRgb segmentation_colors[max_segmentation_colors] = {
  RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
  RNRgb(0.8, 0.5, 0.5), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
  RNRgb(0.4, 0.3, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
  RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),
};



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////


// GLUT variables 

static int GLUTinitialized = 0;
static int GLUTwindow = 0;
static int GLUTwindow_height = 480;
static int GLUTwindow_width = 640;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTmouse_down[2] = { 0, 0 };
static int GLUTmouse_movement = 0;
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// Program variables

static GSVScene *scene = NULL;
static RNArray<StreetSide *> streetsides;
static RNArray<Keystroke *> keystrokes;
static StreetSide *selected_streetside = NULL;
static int selected_streetside_index = -1;
static Channel *selected_channel = NULL;
static int selected_channel_index = 0;
static RNArray<ObjectSegmentation *> selected_segmentations;
static R2Point selected_grid_position(RN_UNKNOWN, RN_UNKNOWN);
static R2Point selection_box_points[2] = { R2Point(0, 0), R2Point(0, 0) };
static RNBoolean selection_box_active = FALSE;
static R2Point split_line_points[2] = { R2Point(0, 0), R2Point(0, 0) };
static RNBoolean split_line_active = FALSE;
static R2Box viewing_window = R2null_box;
static RNInterval viewing_range(0,1);
static RNAngle viewing_angle = 0;
static int color_scheme = RGB_COLOR_SCHEME; 
static int display_type = 0; // 0=channel, 1=points
static int show_channels_in_rgb = 0; 
static int show_segmentation_boundaries = 1; 
static int show_labels = 1; 
static int show_image = 0; 



////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

static RNRgb 
Color(RNScalar value)
{
  // Check for unknown value
  if (value == R2_GRID_UNKNOWN_VALUE) {
    if (show_channels_in_rgb == 0) return RNRgb(1, 0.5, 0);
    else return RNblack_rgb;
  }

  // Normalize color
  RNScalar value_min = viewing_range.Min();
  RNScalar value_width = viewing_range.Diameter();
  RNScalar value_scale = (value_width > 0) ? 1.0 / value_width : 1.0;
  RNScalar normalized_value = value_scale * (value - value_min);

  // Compute color
  RNRgb c(0, 0, 0);
  if (show_channels_in_rgb == 0) {
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



static int
LoadColor(const StreetSide *streetside, const Channel *channel, int i, int j, int color_scheme = RGB_COLOR_SCHEME)
{
  // Useful constant
  RNRgb unknown_color(1.0, 0.5, 0.0);

  // Check color scheme
  if (color_scheme == RGB_COLOR_SCHEME) {
    // Select color based on point colors
    RNRgb color(unknown_color);
    if (streetside && streetside->scanline_grid && streetside->point_backpointers) {
      int grid_index;
      streetside->scanline_grid->IndicesToIndex(i, j, grid_index);
      ObjectPoint *point = streetside->point_backpointers[grid_index];
      if (point) color = point->Color();
    }
    RNLoadRgb(color);
    return 1;
  }
  else if (color_scheme == CHANNEL_COLOR_SCHEME) {
    // Select color based on channel grid values
    if (channel->grids.NEntries() == 1) {
      RNRgb color = Color(channel->grids[0]->GridValue(i, j));
      RNLoadRgb(color);
      return 1;
    }
    else if (channel->grids.NEntries() == 2) {
      RNScalar red = channel->grids[0]->GridValue(i, j);
      if (red != R2_GRID_UNKNOWN_VALUE) {
        RNScalar green = channel->grids[1]->GridValue(i, j);
        if (green != R2_GRID_UNKNOWN_VALUE) {
          RNLoadRgb(red, green, 0.0);
          return 1;
        }
      }
    }
  }
  else if (color_scheme == SIMILARITY_COLOR_SCHEME) {
    // Select color based on similarity to selected grid position
    if (streetside->point_backpointers && streetside->scanline_grid) {
      if ((selected_grid_position.X() != RN_UNKNOWN) && (selected_grid_position.Y() != RN_UNKNOWN)) {
        int selected_i = (int) (selected_grid_position.X() + 0.5);
        if ((selected_i >= 0) && (selected_i < streetside->scanline_grid->XResolution())) {
          int selected_j = (int) (selected_grid_position.Y() + 0.5);
          if ((selected_j >= 0) && (selected_j < streetside->scanline_grid->YResolution())) {
            int selected_grid_index;
            streetside->scanline_grid->IndicesToIndex(selected_i, selected_j, selected_grid_index);
            ObjectPoint *selected_point = streetside->point_backpointers[selected_grid_index];
            if (selected_point) {
              const ObjectDescriptor& selected_descriptor = selected_point->Descriptor();
              if (selected_descriptor.NValues() > 0) {
                int grid_index;
                streetside->scanline_grid->IndicesToIndex(i, j, grid_index);
                ObjectPoint *point = streetside->point_backpointers[grid_index];
                if (point) {
                  const ObjectDescriptor& descriptor = point->Descriptor();
                  if ((descriptor.NValues() > 0) && (descriptor.NValues() == selected_descriptor.NValues())) {
                    double d = descriptor.Distance(selected_descriptor);
                    RNRgb color = Color(d/descriptor.NValues());
                    RNLoadRgb(color);
                    return 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else if (color_scheme == SEGMENTATION_COLOR_SCHEME) {
    // Select color based on segmentation index
    if (streetside && streetside->scanline_grid && streetside->segmentation_backpointers) {
      int grid_index;
      streetside->scanline_grid->IndicesToIndex(i, j, grid_index);
      ObjectSegmentation *segmentation = streetside->segmentation_backpointers[grid_index];
      if (segmentation) {
        int segmentation_index = segmentation->ParseIndex();
        RNLoadRgb(segmentation_colors[segmentation_index % max_segmentation_colors]);
        return 1;
      }
    }
  }
  else if (color_scheme == LABEL_COLOR_SCHEME) {
    // Select color based on label color
    RNRgb color(1.0, 0.5, 0);
    if (streetside && streetside->scanline_grid && streetside->segmentation_backpointers) {
      int grid_index;
      streetside->scanline_grid->IndicesToIndex(i, j, grid_index);
      ObjectSegmentation *segmentation = streetside->segmentation_backpointers[grid_index];
      if (segmentation) {
        ObjectLabel *label = segmentation->GroundTruthLabel();
        if (label) {
          RNLoadRgb(label->Color());
          return 1;
        }
      }
    }
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    // Select color based on grid index
    unsigned char rgba[4];
    int grid_index;
    channel->grids[0]->IndicesToIndex(i, j, grid_index);
    int item_buffer_value = grid_index + 1;
    rgba[0] = item_buffer_value & 0xFF;
    rgba[1] = (item_buffer_value >> 8) & 0xFF;
    rgba[2] = (item_buffer_value >> 16) & 0xFF;
    rgba[3] = (item_buffer_value >> 24) & 0xFF;
    glColor4ubv(rgba);
    return 1;
  }

  // Did not load color
  return 0;
  }



static void 
DrawText(double x, double y, const char *s)
{
  // Draw text string s and position p
  glRasterPos2d(x, y);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}
  


static void 
DrawText(const R2Point& p, const char *s)
{
  // Draw text string s and position p
  DrawText(p[0], p[1], s);
}
  


static void 
DrawGSVImage(GSVImage *image)
{
  // Persistent variables
  static GLuint texture_id = 0;
  static GSVImage *last_image = NULL;

  // Load texture
  if (image != last_image) {
    last_image = image;
    if (texture_id > 0) glDeleteTextures(1, &texture_id);
    texture_id = 0;

    // Create texture
    R2Image *texture = image->UndistortedImage();
    if (texture) {
      // Create texture
      glGenTextures(1, &texture_id);
      glBindTexture(GL_TEXTURE_2D, texture_id);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      assert(texture_id != 0);
      delete texture;
    }
  }

  // Draw textured polygon
  if (texture_id > 0) {
    // Enable texture
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glEnable(GL_TEXTURE_2D);

    // Set projection matrix
    glMatrixMode(GL_PROJECTION);  
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

    // Set model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    double scale = (double) GLUTwindow_height / (double) image->Height();
    scale *= 0.5;
    // glTranslated(GLUTwindow_width/2, GLUTwindow_height/2, 0);
    glTranslated(GLUTwindow_width - scale*image->Width(), GLUTwindow_height - scale*image->Height(), 0);
    glScaled(scale, scale, 1);
    // glTranslated(-image->Width()/2, -image->Height()/2, 0);

    // Draw textured polygon
    glColor3d(1,1,1);
    glBegin(GL_POLYGON);
    glTexCoord2d(0, 0);
    glVertex2d(0.0, 0.0);
    glTexCoord2d(1, 0);
    glVertex2d(image->Width()-1, 0.0);
    glTexCoord2d(1, 1);
    glVertex2d(image->Width()-1, image->Height()-1);
    glTexCoord2d(0, 1);
    glVertex2d(0.0, image->Height()-1);
    glEnd();

    // Reset projection matrix
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    // Reset model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Disable texture
    glDisable(GL_TEXTURE_2D);
  }
}



static void 
DrawSelectionBox(void) 
{
  // Set rendering modes
  glDrawBuffer(GL_FRONT);
  glLineStipple(1, 0xFF00);
  glEnable(GL_LINE_STIPPLE);
  glLogicOp(GL_XOR);
  glEnable(GL_COLOR_LOGIC_OP);
  glLineWidth(3);

  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Draw box
  glBegin(GL_LINE_LOOP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex2f(selection_box_points[0][0], selection_box_points[0][1]);
  glVertex2f(selection_box_points[0][0], selection_box_points[1][1]);
  glVertex2f(selection_box_points[1][0], selection_box_points[1][1]);
  glVertex2f(selection_box_points[1][0], selection_box_points[0][1]);
  glEnd();

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Reset rendering modes
  glLineWidth(1);
  glDisable(GL_LINE_STIPPLE);
  glDisable(GL_COLOR_LOGIC_OP);
  glDrawBuffer(GL_BACK);
  glFlush();
}



static void 
DrawSplitLine(void) 
{
  // Set rendering modes
  glDrawBuffer(GL_FRONT);
  glLineStipple(1, 0xF0F0);
  glEnable(GL_LINE_STIPPLE);
  glLogicOp(GL_XOR);
  glEnable(GL_COLOR_LOGIC_OP);
  glLineWidth(3);

  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Draw line
  glBegin(GL_LINES);
  glColor3d(1.0, 1.0, 1.0);
  glVertex2f(split_line_points[0][0], split_line_points[0][1]);
  glVertex2f(split_line_points[1][0], split_line_points[1][1]);
  glEnd();

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Reset rendering modes
  glLineWidth(1);
  glDisable(GL_LINE_STIPPLE);
  glDisable(GL_COLOR_LOGIC_OP);
  glDrawBuffer(GL_BACK);
  glFlush();
}



static void
DrawScene(void)
{
  // Draw selected streetside
  StreetSide *streetside = streetsides.Kth(selected_streetside_index);
  if (show_labels) streetside->Draw(LABEL_COLOR_SCHEME);
  streetside->Draw(color_scheme);
  if (show_labels) streetside->Draw(LABEL_COLOR_SCHEME);
}



static void
DrawKeystrokes(void)
{
  // Check point type
  if (display_type != 0) return;

  // Compute offset 
  double font_width = 8;
  double font_height = 10;
  double dx = -0.5 * font_width * viewing_window.XLength() / (double) GLUTwindow_width;
  double dy = -0.5 * font_height * viewing_window.YLength() / (double) GLUTwindow_height;
  R2Vector offset(dx, dy);

  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(viewing_window.XMin(), viewing_window.XMax(), viewing_window.YMin(), viewing_window.YMax()); 

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Draw keystrokes
  RNLoadRgb(RNyellow_rgb);
  for (int i = 0; i < keystrokes.NEntries(); i++) {
    Keystroke *keystroke = keystrokes.Kth(i);
    static char buffer[2] = { 0, 0 };
    buffer[0] = keystroke->key;
    DrawText(keystroke->p + offset, buffer);
  }

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}



////////////////////////////////////////////////////////////////////////
// Pick functions
////////////////////////////////////////////////////////////////////////

static int 
Pick(int x, int y, 
  R2Point *picked_grid_position = NULL, R3Point *picked_world_position = NULL, 
  StreetSide **picked_streetside = NULL)
{
  // Get/check useful variables
  StreetSide *streetside = selected_streetside;
  if (!streetside) return 0;
  const R2Grid *grid = streetside->scanline_grid;

  // Get grid index at (x,y)
  int grid_index = -1;
  if (display_type == 0) {
    int ix = (int) (x * viewing_window.XLength() / (double) GLUTwindow_width + viewing_window.XMin() + 0.5);
    int iy = (int) (y * viewing_window.YLength() / (double) GLUTwindow_height + viewing_window.YMin() + 0.5);
    grid->IndicesToIndex(ix, iy, grid_index);
  }
  else {        
    // Draw scene with pick colors
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    streetside->Draw(PICK_COLOR_SCHEME);
    glFinish();

    // Read color buffer at cursor position
    unsigned char rgba[4];
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
    unsigned int r = rgba[0] & 0xFF;
    unsigned int g = rgba[1] & 0xFF;
    unsigned int b = rgba[2] & 0xFF;
    unsigned int a = rgba[3] & 0xFF;
    int item_buffer_value = (a << 24) | (b << 16) | (g << 8) | r;
    if (item_buffer_value == 0) return 0;

    // Determine/check grid index
    grid_index = item_buffer_value - 1;
    if (grid_index < 0) return 0;
    if (grid_index >= grid->NEntries()) return 0;
  }

  // Return streetside
  if (picked_streetside) {
    *picked_streetside = streetside;
  }

  // Return grid position
  if (picked_grid_position) {
    int ix, iy;
    grid->IndexToIndices(grid_index, ix, iy);
    R2Point grid_position(ix, iy);
    *picked_grid_position = grid_position;
  }

  // Return position
  if (picked_world_position) {
#if 1
    // Find hit position
    ObjectPoint *point = streetside->point_backpointers[grid_index];
    if (point) *picked_world_position = point->Position();
    else return 0;
#else
    // Find hit position
    GLfloat depth;
    GLdouble p[3];
    GLint viewport[4];
    GLdouble modelview_matrix[16];
    GLdouble projection_matrix[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    gluUnProject(x, y, depth, modelview_matrix, projection_matrix, viewport, &(p[0]), &(p[1]), &(p[2]));
    R3Point world_position(p[0], p[1], p[2]);
    *picked_world_position = world_position;
#endif
  }
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Selection functions
////////////////////////////////////////////////////////////////////////

static int 
EmptySelections(void)
{
  selected_segmentations.Empty();
  return 1;
}



static int
UpdateSelections(ObjectSegmentation *segmentation, int operation)
{
  // Update selection
  if (operation == 0) {
    // Replace
    EmptySelections();
    selected_segmentations.Insert(segmentation);
  }
  else if (operation == 1) {
    // Add 
    RNArrayEntry *entry = selected_segmentations.FindEntry(segmentation);
   if (!entry) selected_segmentations.Insert(segmentation);
  }
  else if (operation == 2) {
    // Remove 
    RNArrayEntry *entry = selected_segmentations.FindEntry(segmentation);
    if (entry) selected_segmentations.RemoveEntry(entry);
  }
  else {
    // Toggle
    RNArrayEntry *entry = selected_segmentations.FindEntry(segmentation);
    if (entry) selected_segmentations.RemoveEntry(entry);
    else selected_segmentations.Insert(segmentation);
  }

  // Return success
  return 1;
}



static int
UpdateSelections(int cx, int cy, int operation)
{
  // Empty selections 
  if (operation == 0) { EmptySelections(); operation = 1; }

  // Pick grid position
  R2Point p;
  if (!Pick(cx, cy, &p)) return 0;

  // Get segmentation at ix, iy
  int grid_index;
  int ix = (int) (p.X() + 0.5);
  int iy = (int) (p.Y() + 0.5);
  if (!selected_streetside || !selected_streetside->scanline_grid) return 0;
  selected_streetside->scanline_grid->IndicesToIndex(ix, iy, grid_index);
  ObjectSegmentation *segmentation = selected_streetside->segmentation_backpointers[grid_index];
  if (!segmentation) return 0;

  // Update selections
  return UpdateSelections(segmentation, operation);
}



static int
UpdateSelections(const R2Box& selection_box, int operation = 0)
{
  // Empty selections 
  if (operation == 0) { EmptySelections(); operation = 1; }

  // Check stuff
  if (!selected_streetside) return 0;
  R2Grid *scanline_grid = selected_streetside->scanline_grid;
  if (!scanline_grid) return 0;
  ObjectSegmentation **segmentation_backpointers = selected_streetside->segmentation_backpointers;
  if (!segmentation_backpointers) return 0;
  ObjectParse *parse = selected_streetside->parse;
  if (!parse) return 0;
  R3Viewer viewer = selected_streetside->PointViewer();

  // Draw scene with pick colors
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  selected_streetside->Draw(PICK_COLOR_SCHEME);
  glFinish();

  // Determine extent of selection box
  int xmin = (int) (selection_box.XMin() + 0.5);
  int ymin = (int) (selection_box.YMin() + 0.5);
  int width = (int) (selection_box.XLength());
  int height = (int) (selection_box.YLength());
  if (xmin < 0) xmin = 0;
  if (xmin > GLUTwindow_width - 1) xmin = GLUTwindow_width - 1;
  if (ymin < 0) ymin = 0;
  if (ymin > GLUTwindow_height - 1) ymin = GLUTwindow_height - 1;
  if (width > GLUTwindow_width - xmin - 1) width = GLUTwindow_width - xmin - 1;
  if (height > GLUTwindow_height - ymin - 1) height = GLUTwindow_height - ymin - 1;
  if ((width <= 0) || (height <= 0)) return 0;

  // Read color buffer within selection box
  unsigned char *pixels = new unsigned char [4 * width * height ];
  glReadPixels(xmin, ymin, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

  // Allocate marks
  int *marks = new int [ parse->NDetections() ];
  for (int i = 0; i < parse->NDetections(); i++) marks[i] = 0;

  // Consider segmentations within selection box
  for (int i = 0; i < width * height; i++) {
    unsigned char *rgba = pixels + i*4;
    unsigned int r = rgba[0] & 0xFF;
    unsigned int g = rgba[1] & 0xFF;
    unsigned int b = rgba[2] & 0xFF;
    unsigned int a = rgba[3] & 0xFF;
    int item_buffer_value = (a << 24) | (b << 16) | (g << 8) | r;
    if (item_buffer_value == 0) continue;

    // Get segmentation
    int grid_index = item_buffer_value - 1;
    if (grid_index < 0) continue;
    if (grid_index >= scanline_grid->NEntries()) continue;
    ObjectSegmentation *segmentation = segmentation_backpointers[grid_index];
    if (!segmentation) continue;

    // Check if segmentation has already been considered
    ObjectDetection *detection = segmentation->Detection();
    if (!detection) continue;
    int detection_index = detection->ParseIndex();
    if (marks[detection_index]) continue;
    marks[detection_index] = 1;

    // Check if segmentation is entirely within selection box
    RNBoolean inside = TRUE;
    const ObjectPointSet& pointset = segmentation->PointSet();
    for (int i = 0; i < pointset.NPoints(); i++) {
      ObjectPoint *point = pointset.Point(i);
      if (display_type == 0) {
        int ix, iy;
        int grid_index = (int) (point->Value() + 0.5);
        selected_streetside->scanline_grid->IndexToIndices(grid_index, ix, iy);
        double sx = (ix - viewing_window.XMin()) * GLUTwindow_width / viewing_window.XLength();
        double sy = (iy - viewing_window.YMin()) * GLUTwindow_height / viewing_window.YLength();
        R2Point screen_position(sx, sy);
        if (!R2Contains(selection_box, screen_position)) {
          inside = FALSE;
          break;
        }
      }
      else {
        R2Point screen_position = viewer.ViewportPoint(point->Position());
        if (!R2Contains(selection_box, screen_position)) {
          inside = FALSE;
          break;
        }
      }
    }

    // Update selections
    if (inside) UpdateSelections(segmentation, operation);
  }

  // Delete temporary data
  delete [] pixels;
  delete [] marks;

  // Return success
  return 1;
}


////////////////////////////////////////////////////////////////////////
// Streetside and Channel selection functions
////////////////////////////////////////////////////////////////////////

static void
SelectChannel(int channel_index)
{
  // Check selected streetside
  if (!selected_streetside) return;
  assert(!selected_streetside->channels.IsEmpty());

  // Check index
  if (channel_index < 0) channel_index = 0;
  if (channel_index > selected_streetside->channels.NEntries()-1) channel_index =  selected_streetside->channels.NEntries()-1;
  Channel *channel = selected_streetside->channels.Kth(channel_index);
  assert(channel->streetside == selected_streetside);

  // Update display variables
  if (!selected_channel ||
      (selected_channel->width != channel->width) ||
      (selected_channel->height != channel->height)) {
    RNScalar window_aspect = (double) GLUTwindow_width / (double) GLUTwindow_height;
    RNScalar channel_aspect = (double) channel->width / (double) channel->height;
    R2Point origin(channel->width/2, channel->height/2);
    R2Vector diagonal(channel->width/2, channel->height/2);
    diagonal[0] *= window_aspect / channel_aspect;
    viewing_window = R2Box(origin - diagonal, origin + diagonal);
  }

  // Update min and max values
  viewing_range = channel->grids[0]->Range();

  // Update selected channel 
  selected_channel_index = channel_index;
  selected_channel = channel;

  // Set window title
  if (GLUTinitialized) {
    char title[4096];
    GSVScan *scan = selected_streetside->scan;
    int scan_index = scan->SegmentIndex();
    GSVSegment *segment = scan->Segment();
    int segment_index = segment->RunIndex();
    GSVRun *run = segment->Run();
    sprintf(title, "%s %02d %02d %s %s", run->Name(), segment_index, scan_index, channel->parameterization, channel->name);
    glutSetWindowTitle(title);
  }
}



static void 
SelectStreetSide(int streetside_index)
{
  // Check index
  if (streetside_index < 0) streetside_index = 0;
  if (streetside_index > streetsides.NEntries()-1) streetside_index =  streetsides.NEntries()-1;

  // Check if no change
  if (streetside_index == selected_streetside_index) return;

  // Empty selections
  EmptySelections();

  // Unselect old streetside
  if (selected_streetside) {
    // Release stuff for old streetside
    selected_streetside->Release();

    // Unselect old streetside
    selected_streetside_index = -1;
    selected_streetside = NULL;

    // Unselect channel
    selected_channel = NULL;
  }

  // Select new streetside
  if (streetside_index >= 0) {
    // Select streetside
    selected_streetside_index = streetside_index;
    selected_streetside = streetsides.Kth(streetside_index);

    // Read stuff for selected streetside
    if (!selected_streetside->Read()) {
      fprintf(stderr, "Unable to read streetside\n");
      exit(-1);
    }

    // Select channel
    if (selected_channel_index >= 0) {
      SelectChannel(selected_channel_index);
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Segementation functions
////////////////////////////////////////////////////////////////////////

static int
MergeSegmentation(ObjectSegmentation *segmentation0, ObjectSegmentation *segmentation1)
{
  // Move points from segmentation1 into segmentation0
  const ObjectPointSet& pointset1 = segmentation1->PointSet();
  for (int i = 0; i < pointset1.NPoints(); i++) {
    ObjectPoint *point1 = pointset1.Point(i);
    segmentation0->pointset.InsertPoint(*point1);
  }

  // Update backpointers 
  for (int i = 0; i < segmentation0->PointSet().NPoints(); i++) {
    ObjectPoint *point = segmentation0->PointSet().Point(i);
    int grid_index = (int) (point->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation0;
    selected_streetside->point_backpointers[grid_index] = point;
  }

  // Delete segmentation1
  delete segmentation1;

  // Return success
  return 1;
}



static int
MergeSelectedSegmentations(void)
{
  // Get first selected segmentation
  ObjectSegmentation *segmentation = selected_segmentations.Kth(0);

  // Merge other selected segmentations
  for (int i = 1; i < selected_segmentations.NEntries(); i++) {
    MergeSegmentation(segmentation, selected_segmentations[i]);
  }

  // Reset selections
  UpdateSelections(segmentation, 0);

  // Return success
  return 1;
}



#if 0
static int
SplitSegmentation(ObjectSegmentation *segmentation, 
 const RNArray<ObjectPoint *>& outlier_points,
 ObjectSegmentation **outlier_segmentation = NULL)
{
  // Check stuff
  if (!segmentation) return 0;
  if (outlier_points.NEntries() == 0) return 0;
  const RNScalar mark_value = -333;

  // Allocate temporary memory to save point values
  RNScalar *saved_values = new RNScalar [ outlier_points.NEntries() ];

  // Mark point values in pointset2
  for (int i = 0; i < outlier_points.NEntries(); i++) {
    ObjectPoint *point = outlier_points.Kth(i);
    saved_values[i] = point->Value();
    assert(point->Value() != mark_value);
    point->SetValue(mark_value);
  }

  // Create pointset1 
  ObjectPointSet pointset1;
  const ObjectPointSet& pointset = segmentation->PointSet();
  for (int i = 0; i < pointset.NPoints(); i++) {
    ObjectPoint *point = pointset.Point(i);
    if (point->Value() == mark_value) continue;
    pointset1.InsertPoint(*point);
  }

  // Check pointset1
  if (pointset1.NPoints() == 0) return 0;
  if (pointset1.NPoints() == pointset.NPoints()) return 0;

  // Create pointset2
  ObjectPointSet pointset2;
  for (int i = 0; i < outlier_points.NEntries(); i++) {
    ObjectPoint *point = outlier_points.Kth(i);
    point->SetValue(saved_values[i]);
    pointset2.InsertPoint(*point);
  }

  // Check pointset2
  if (pointset2.NPoints() == 0) return 0;
  if (pointset2.NPoints() == pointset.NPoints()) return 0;
  int swap = (pointset2.NPoints() > pointset1.NPoints()) ? 1 : 0;

  // Update segmentation1
  ObjectSegmentation *segmentation1 = segmentation;
  if (swap) segmentation1->SetPointSet(pointset2);
  else segmentation1->SetPointSet(pointset1);

  // Create segmentation2
  ObjectDetection *detection2 = new ObjectDetection();
  selected_streetside->parse->InsertDetection(detection2);
  ObjectSegmentation *segmentation2 = new ObjectSegmentation();
  detection2->InsertSegmentation(segmentation2);
  if (swap) segmentation2->SetPointSet(pointset1);
  else segmentation2->SetPointSet(pointset2);

  // Update backpointers for segmentation1
  for (int i = 0; i < segmentation1->PointSet().NPoints(); i++) {
    ObjectPoint *point = segmentation1->PointSet().Point(i);
    int grid_index = (int) (point->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation1;
    selected_streetside->point_backpointers[grid_index] = point;
  }

  // Update backpointers for segmentation1
  for (int i = 0; i < segmentation2->PointSet().NPoints(); i++) {
    ObjectPoint *point = segmentation2->PointSet().Point(i);
    int grid_index = (int) (point->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation2;
    selected_streetside->point_backpointers[grid_index] = point;
  }

  // Delete temporary memory
  delete [] saved_values;

  // Return new segmentation
  if (outlier_segmentation) *outlier_segmentation = segmentation2;

  // Return success
  return 1;
}



static int
SplitSegmentation(ObjectSegmentation *segmentation, 
 const R3Plane& plane, ObjectSegmentation **outlier_segmentation = NULL)
{
  // Check stuff
  if (!segmentation) return 0;

  // Create pointsets
  ObjectPointSet pointset1, pointset2;
  const ObjectPointSet& pointset = segmentation->PointSet();
  for (int i = 0; i < pointset.NPoints(); i++) {
    ObjectPoint *point = pointset.Point(i);
    const R3Point& position = point->Position();
    if (R3SignedDistance(plane, position) < 0) pointset2.InsertPoint(*point);
    else pointset1.InsertPoint(*point);
  }

  // Check pointsets
  if (pointset1.NPoints() == 0) return 0;
  if (pointset1.NPoints() == pointset.NPoints()) return 0;
  if (pointset2.NPoints() == 0) return 0;
  if (pointset2.NPoints() == pointset.NPoints()) return 0;

  // Update segmentation1
  ObjectSegmentation *segmentation1 = segmentation;
  segmentation1->SetPointSet(pointset1);

  // Create segmentation2
  ObjectDetection *detection2 = new ObjectDetection();
  selected_streetside->parse->InsertDetection(detection2);
  ObjectSegmentation *segmentation2 = new ObjectSegmentation();
  detection2->InsertSegmentation(segmentation2);
  segmentation2->SetPointSet(pointset2);

  // Update backpointers
  for (int i = 0; i < pointset2.NPoints(); i++) {
    ObjectPoint *point2 = pointset2.Point(i);
    int grid_index = (int) (point2->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation2;
  }

  // Return new segmentation
  if (outlier_segmentation) *outlier_segmentation = segmentation2;

  // Return success
  return 1;
}
#endif



static int
SplitSegmentation(ObjectSegmentation *segmentation, 
  const R2Line& split_line,
  ObjectSegmentation **outlier_segmentation = NULL)
{
  // Check stuff
  if (!segmentation) return 0;
  if (!selected_streetside) return 0;
  if (!selected_streetside->scanline_grid) return 0;
  R3Viewer viewer = selected_streetside->PointViewer();

  // Create pointsets
  ObjectPointSet pointset1, pointset2;
  const ObjectPointSet& pointset = segmentation->PointSet();
  for (int i = 0; i < pointset.NPoints(); i++) {
    ObjectPoint *point = pointset.Point(i);

    // Get screen position
    R2Point screen_position(0, 0);
    if (display_type == 0) {
      int ix, iy;
      int grid_index = (int) (point->Value() + 0.5);
      selected_streetside->scanline_grid->IndexToIndices(grid_index, ix, iy);
      double sx = (ix - viewing_window.XMin()) * GLUTwindow_width / viewing_window.XLength();
      double sy = (iy - viewing_window.YMin()) * GLUTwindow_height / viewing_window.YLength();
      screen_position.Reset(sx, sy);
    }
    else {
      screen_position = viewer.ViewportPoint(point->Position());
    }

    // Assign point to pointset
    RNScalar d = R2SignedDistance(split_line, screen_position);
    if (d < 0) pointset2.InsertPoint(*point);
    else pointset1.InsertPoint(*point);
  }

  // Check pointsets
  if (pointset1.NPoints() == 0) return 0;
  if (pointset1.NPoints() == pointset.NPoints()) return 0;
  if (pointset2.NPoints() == 0) return 0;
  if (pointset2.NPoints() == pointset.NPoints()) return 0;

  // Update segmentation1
  ObjectSegmentation *segmentation1 = segmentation;
  segmentation1->SetPointSet(pointset1);

  // Create segmentation2
  ObjectDetection *detection2 = new ObjectDetection();
  selected_streetside->parse->InsertDetection(detection2);
  ObjectSegmentation *segmentation2 = new ObjectSegmentation();
  detection2->InsertSegmentation(segmentation2);
  segmentation2->SetPointSet(pointset2);

  // Update backpointers for segmentation1
  for (int i = 0; i < segmentation1->PointSet().NPoints(); i++) {
    ObjectPoint *point = segmentation1->PointSet().Point(i);
    int grid_index = (int) (point->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation1;
    selected_streetside->point_backpointers[grid_index] = point;
  }

  // Update backpointers for segmentation1
  for (int i = 0; i < segmentation2->PointSet().NPoints(); i++) {
    ObjectPoint *point = segmentation2->PointSet().Point(i);
    int grid_index = (int) (point->Value() + 0.5);
    selected_streetside->segmentation_backpointers[grid_index] = segmentation2;
    selected_streetside->point_backpointers[grid_index] = point;
  }

  // Return new segmentation
  if (outlier_segmentation) *outlier_segmentation = segmentation2;

  // Return success
  return 1;
}



static int
SplitSelectedSegmentations(const R2Line& split_line)
{
  // Split each selected segmentation
  for (int i = 0; i < selected_segmentations.NEntries(); i++) {
    ObjectSegmentation *segmentation = selected_segmentations.Kth(i);
    SplitSegmentation(segmentation, split_line);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Assignment functions
////////////////////////////////////////////////////////////////////////

static int
InsertAssignment(int key)
{
  printf("%d\n", key);
  // Check selected segmentations
  if (selected_segmentations.IsEmpty()) return 0;
  printf("1237\n");
  // Get model
  ObjectParse *parse = selected_streetside->parse;
  printf("%p\n", parse);
  ObjectLabel *label = parse->LabelFromKey(key);
  printf("%p\n", label);
  if (!label) return 0;
  printf("%d\n", label->NModels());
  if (label->NModels() == 0) return 0;
  ObjectModel *model = label->Model(0);
  printf("%p\n", model);
  if (!model) return 0;
  printf("1245\n");
  // Merge selected segmentations
  ObjectSegmentation *segmentation = selected_segmentations.Kth(0);
  for (int i = 1; i < selected_segmentations.NEntries(); i++) {
    MergeSegmentation(segmentation, selected_segmentations[i]);
  }
  printf("1251\n");
  // Remove previous assignments
  segmentation->EmptyAssignments();
  printf("1254\n");
  // Create assignment
  ObjectAssignment *assignment = new ObjectAssignment(segmentation, model, R3identity_affine, 1.0, TRUE);
  segmentation->InsertAssignment(assignment);
  model->InsertAssignment(assignment);
  printf("1259\n");
  // Empty selections
  EmptySelections();
  printf("1262\n");
  // Return success
  return 1;
}



static int
RemoveAssignment(void)
{
  // Check selected segmentations
  if (selected_segmentations.IsEmpty()) return 0;

  // Delete assignments on selected segmentations
  for (int i = 0; i < selected_segmentations.NEntries(); i++) {
    ObjectSegmentation *segmentation = selected_segmentations.Kth(i);
    while (segmentation->NAssignments() > 0) {
      ObjectAssignment *assignment = segmentation->Assignment(0);
      delete assignment;
    }
  }

  // Empty selections
  EmptySelections();

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Keystroke functions
////////////////////////////////////////////////////////////////////////

#if 0
static void
InsertKeystroke(int c, int x, int y)
{
  // Pick point
  R2Point p;
  if (Pick(x, y, &p)) {
    // Create and insert keystroke
    Keystroke *keystroke = new Keystroke(selected_streetside, parameterization, p, c);
    keystrokes.Insert(keystroke);
  }
}



static int
FindKeystroke(int c, int x, int y)
{
  // Pick point
  R2Point p;
  int closest_index = -1;
  if (Pick(x, y, &p)) {
    // Search for matching keystroke
    RNLength closest_distance = 2;
    for (int i = 0; i < keystrokes.NEntries(); i++) {
      Keystroke *keystroke = keystrokes.Kth(i);
      if ((c > 0) && (keystroke->key != c)) continue;
      RNLength distance = R2Distance(keystroke->p, p);
      if (distance < closest_distance) {
        closest_distance = distance;
        closest_index = i;
      }
    }
  }

  // Return closest index
  return closest_index;
}



static void
RemoveKeystroke(int c, int x, int y)
{
  // Find keystroke
  int k = FindKeystroke(c, x, y);
  if (k < 0) return;

  // Remove keystroke
  Keystroke *keystroke = keystrokes.Kth(k);
  keystrokes.RemoveKth(k);
  delete keystroke;
}
#endif



////////////////////////////////////////////////////////////////////////
// GSV input/output functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
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
  if (!scene->ReadFile(filename, FALSE)) {
    delete scene;
    return NULL;
  }

  // Create streetsides
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        StreetSide *streetside = new StreetSide(scan);
        streetside->index = streetsides.NEntries();
        assert(streetside->index == scan->SceneIndex());
        streetsides.Insert(streetside);
      }
    }
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
    printf("  # Scans = %d\n", scene->NScans());
    printf("  # Scanlines = %d\n", scene->NScanlines());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



////////////////////////////////////////////////////////////////////////
// Keystroke input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadKeystrokes(GSVScene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open keystroke file %s\n", filename);
    return 0;
  }

  // Read keystrokes
  int count = 0;
  char buffer[4096];
  while (fgets(buffer, 4096, fp)) {
    count++;
    double x, y;
    int segment_index, scan_index, key, dummy; 
    char cmd[1024], run_name[1024], parameterization[1024];
    if ((sscanf(buffer, "%s", cmd) == (unsigned int) 1) && (cmd[0] != '#')) {
      if (sscanf(buffer, "%d%s%d%d%s%d%lf%lf%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", 
        &dummy, run_name, &segment_index, &scan_index, 
        parameterization, &key, &x, &y, 
        &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,
        &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 24) {
        fprintf(stderr, "Unable to read line %d of %s\n", count, filename);
        return 0; 
      }
      GSVRun *run = scene->Run(run_name);
      if (!run) { 
        fprintf(stderr, "Invalid run %s at line %d\n", run_name, count);
        return 0; 
      }
      if ((segment_index < 0) || (segment_index >= run->NSegments())) { 
        fprintf(stderr, "Invalid segment index %d at line %d of %s\n", segment_index, count, filename); 
        return 0; 
      }
      GSVSegment *segment = run->Segment(segment_index);
      if ((scan_index < 0) || (scan_index >= segment->NScans())) { 
        fprintf(stderr, "Invalid scan index %d at line %d of %s\n", scan_index, count, filename); 
        return 0; 
      }
      GSVScan *scan = segment->Scan(scan_index);
      if ((scan->SceneIndex() < 0) || (scan->SceneIndex() >= streetsides.NEntries())) { 
        fprintf(stderr, "Invalid streetside %d at line %d of %s\n", scan->SceneIndex(), count, filename); 
        return 0; 
      }
      StreetSide *streetside = streetsides.Kth(scan->SceneIndex());
      if (!streetside) { 
        fprintf(stderr, "Unable to find streetside at line %d of %s\n", count, filename); 
        return 0; 
      }
      Keystroke *keystroke = new Keystroke(streetside, parameterization, R2Point(x, y), key);
      keystrokes.Insert(keystroke);
    }
  }

  // Close file
  fclose(fp);

  // Print message
  if (print_verbose) {
    printf("Read keystrokes from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Keystrokes = %d\n", keystrokes.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteKeystrokes(GSVScene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open keystroke file %s\n", filename);
    return 0;
  }

  // Write keystrokes
  for (int i = 0; i < keystrokes.NEntries(); i++) {
    Keystroke *keystroke = keystrokes.Kth(i);
    StreetSide *streetside = keystroke->streetside;  
    GSVScan *scan = streetside->scan;
    GSVSegment *segment = scan->Segment();
    GSVRun *run = segment->Run();
    fprintf(fp, "0   %-20s %02d %02d  %4s %4d  %9.3f %9.3f  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n", 
      run->Name(), segment->RunIndex(), scan->SegmentIndex(),
      keystroke->parameterization, keystroke->key, keystroke->p.X(), keystroke->p.Y());
  }

  // Close file
  fclose(fp);

  // Print message
  if (print_verbose) {
    printf("Wrote keystrokes to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Keystrokes = %d\n", keystrokes.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Rectangle I/O functions
////////////////////////////////////////////////////////////////////////

static int
WriteScanRectangles(GSVScene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_selected_streetside_index = selected_streetside_index;
  int count = 0;

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open scan rectangles file %s\n", filename);
    return 0;
  }

  // Write rectangles
  for (int i = 0; i < streetsides.NEntries(); i++) {
    StreetSide *streetside = streetsides.Kth(i);

    // Select streetside
    SelectStreetSide(i);

    // Get useful stuff
    ObjectParse *parse = streetside->parse;
    if (!parse) continue;
    R2Grid *grid = streetside->scanline_grid;
    if (!grid) continue;
    GSVScan *scan = streetside->scan;
    if (!scan) continue;
    int scan_index = scan->SegmentIndex();
    GSVSegment *segment = scan->Segment();
    int segment_index = segment->RunIndex();
    GSVRun *run = segment->Run();

    // Write rectangle for every labeled detection
    for (int j = 0; j < parse->NDetections(); j++) {
      ObjectDetection *detection = parse->Detection(j);
      if (detection->NSegmentations() == 0) continue;
      ObjectSegmentation *segmentation = detection->Segmentation(0);
      if (!segmentation) continue;
      ObjectAssignment *assignment = segmentation->GroundTruthAssignment();
      if (!assignment) assignment = segmentation->BestAssignment();
      if (!assignment) continue;
      ObjectModel *model = assignment->Model();
      if (!model) continue;
      ObjectLabel *label = model->Label();
      if (!label) continue;
      const ObjectPointSet& pointset = segmentation->PointSet();
      if (pointset.NPoints() == 0) continue;

      // Compute box
      R2Box box = R2null_box;
      for (int k = 0; k < pointset.NPoints(); k++) {
        int ix, iy;
        ObjectPoint *point = pointset.Point(k);
        int grid_index = (int) (point->Value() + 0.5);
        if ((grid_index < 0) || (grid_index >= grid->NEntries())) continue;
        grid->IndexToIndices(grid_index, ix, iy);
        box.Union(R2Point(ix, iy));
      }

      // Write box
      fprintf(fp, "R  %20s %02d %02d %2s   %s %d %d   %d %d   %g %d   %g %g %g %g\n", 
        run->Name(), segment_index, scan_index, parameterization,
        label->Name(), label->ParseIndex(), model->LabelIndex(), 
        detection->ParseIndex(), segmentation->DetectionIndex(),
        assignment->Score(), (assignment->IsGroundTruth()) ? 1 : 0,
        box.XMin(), box.YMin(), box.XMax(), box.YMax());
    }
  }

  // Close file
  fclose(fp);

  // Restore selected streetside
  if (saved_selected_streetside_index >= 0) {
    SelectStreetSide(saved_selected_streetside_index);
  }

  // Print message
  if (print_verbose) {
    printf("Wrote scan rectangles to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Rectangles = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteImageRectangles(GSVScene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open scan rectangles file %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Print message
  if (print_verbose) {
    printf("Wrote scan rectangles to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Rectangles = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// GLUT callback functions
////////////////////////////////////////////////////////////////////////

void AtExit(void)
{
  // Write keystroke file
  if (output_keystroke_name) {
    WriteKeystrokes(scene, output_keystroke_name);
  }

  // Write scan rectangle file
  if (output_scan_rectangle_name) {
    WriteScanRectangles(scene, output_scan_rectangle_name);
  }

  // Write rectangle file
  if (output_image_rectangle_name) {
    WriteImageRectangles(scene, output_image_rectangle_name);
  }

  // Release streetside (writes files)
  if (selected_streetside) {
    selected_streetside->Release();
  }
}


void GLUTStop(void)
{
  // Exit
  exit(0);
}



void GLUTRedraw(void)
{
  // Check selected streetside
  if (!scene) return;
  if (!selected_streetside) return;

  // Clear window 
  glClearColor(1.0, 0.5, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw scene
  DrawScene();

  // Draw keystrokes
  DrawKeystrokes();

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

  // Update viewing window
  if (selected_streetside) {
    if (selected_channel) {    
      RNScalar window_aspect = (double) GLUTwindow_width / (double) GLUTwindow_height;
      RNScalar channel_aspect = (double) selected_channel->width / (double) selected_channel->height;
      R2Point origin(selected_channel->width/2, selected_channel->height/2);
      R2Vector diagonal(selected_channel->width/2, selected_channel->height/2);
      diagonal[0] *= window_aspect / channel_aspect;
      viewing_window = R2Box(origin - diagonal, origin + diagonal);
    }
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
  
  // Check for drag
  RNBoolean drag = (GLUTmouse_movement > 10);

  // View manipulation
  if (GLUTbutton[0]) {
    if (drag) {
      if (glutGetModifiers() & GLUT_ACTIVE_ALT) {
        if (split_line_active) {
          // Update split line
          DrawSplitLine();
          split_line_points[1] = R2Point(x, y);
          DrawSplitLine();
        }
        else {
          split_line_active = TRUE;
          split_line_points[0] = R2Point(GLUTmouse_down[0], GLUTmouse_down[1]);
          split_line_points[1] = R2Point(x, y);
          DrawSplitLine();
        }
      }
      else {
        if (selection_box_active) {
          // Update selection box
          DrawSelectionBox();
          selection_box_points[1] = R2Point(x, y);
          DrawSelectionBox();
        }
        else {
          selection_box_active = TRUE;
          selection_box_points[0] = R2Point(GLUTmouse_down[0], GLUTmouse_down[1]);
          selection_box_points[1] = R2Point(x, y);
          DrawSelectionBox();
        }
      }
    }
  }
  else if (GLUTbutton[1]) {
    // Zoom and rotate
    RNScalar scale_factor = 1;
    scale_factor *= 1.0-(double)dx/(double)GLUTwindow_width;
    scale_factor *= scale_factor;
    viewing_window.Inflate(scale_factor);
    viewing_angle += 0.01 * dy;
    if (viewing_angle > 0) viewing_angle = 0;
    if (viewing_angle < -0.5*RN_PI) viewing_angle = -0.5*RN_PI;
    glutPostRedisplay();
  }
  else if (GLUTbutton[2]) {
    // Pan and rotate
    RNScalar tx = -dx * viewing_window.XLength() / (double) GLUTwindow_width;
    RNScalar ty = -dy * viewing_window.YLength() / (double) GLUTwindow_height;
    viewing_window.Translate(R2Vector(tx, ty));
    glutPostRedisplay();
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Update mouse movement
  GLUTmouse_movement += dx*dx + dy*dy;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Check if up or down
  if (state == GLUT_UP) {
    // Check for drag
    RNBoolean drag = (GLUTmouse_movement > 10);

    // Check for double click
    static RNBoolean double_click = FALSE;
    static RNTime last_mouse_down_time;
    double_click = !drag && !double_click && (last_mouse_down_time.Elapsed() < 0.4);
    last_mouse_down_time.Read();

    // Process mouse button event
    if (button == 0) {
      if (selection_box_active) {
        // Select segmentations in box
        if (!R2Contains(selection_box_points[0], selection_box_points[1])) {
          R2Box selection_box = R2null_box;
          selection_box.Union(selection_box_points[0]);
          selection_box.Union(selection_box_points[1]);
          int operation = 0;
          if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) operation = 1;
          if (glutGetModifiers() & GLUT_ACTIVE_CTRL) operation = 2;
          UpdateSelections(selection_box, operation);
        }
        DrawSelectionBox();
        selection_box_points[0] = R2zero_point;
        selection_box_points[1] = R2zero_point;
        selection_box_active = FALSE;
        glutPostRedisplay();
      }
      else if (split_line_active) {
        // Split segmentations
        if (!R2Contains(split_line_points[0], split_line_points[1])) {
          R2Line split_line(split_line_points[0], split_line_points[1]);
          SplitSelectedSegmentations(split_line);
        }
        DrawSplitLine();
        split_line_points[0] = R2zero_point;
        split_line_points[1] = R2zero_point;
        split_line_active = FALSE;
        glutPostRedisplay();
      }
      else {
        // Update selections
        int operation = 0;
        if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) operation = 1;
        if (glutGetModifiers() & GLUT_ACTIVE_CTRL) operation = 2;
        UpdateSelections(x, y, operation);
        glutPostRedisplay();
      }
    }
    else if (button == 1) {
      if (!drag) {
        // Query
        Pick(x, y, &selected_grid_position);
        glutPostRedisplay();
      }
    }
    else if ((button == 3) || (button == 4)) {
      // Zoom with wheel
      RNScalar scale_factor = (button == 3) ? 0.9 : 1.1;
      viewing_window.Inflate(scale_factor);
      glutPostRedisplay();
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Reset mouse movement
  GLUTmouse_movement = 0;

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember mouse down position 
  if (state == GLUT_DOWN) {
    GLUTmouse_down[0] = x;
    GLUTmouse_down[1] = y;
  }

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
  case GLUT_KEY_PAGE_DOWN: {
    printf("Changing streetside:\n");
    // Select new streetside
    int d = (key == GLUT_KEY_PAGE_UP) ? -1 : 1;
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      SelectStreetSide(selected_streetside_index-d);
    }
    else {
      SelectChannel(selected_channel_index-d);
    }
    glutPostRedisplay();
    break; }

  case GLUT_KEY_LEFT: 
  case GLUT_KEY_RIGHT: { 
    // Pan view
    int dx = (key == GLUT_KEY_LEFT) ? -20 : 20;
    viewing_window.Translate(R2Vector(dx, 0));
    glutPostRedisplay();
    break; }

  case GLUT_KEY_UP: 
  case GLUT_KEY_DOWN: {
    // Zoom view
    RNScalar scale_factor = (key == GLUT_KEY_UP) ? 0.9 : 1.1;
    viewing_window.Inflate(scale_factor);
    glutPostRedisplay();
    break; }
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
  if ((key >= 'A') && (key <= 'z')) {
    // InsertKeystroke(tolower(key), x, y);
    InsertAssignment(key);
  }
  else {
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
      SelectChannel(key - '1');
      break;
      
    case 2: // ctrl-B
      show_segmentation_boundaries = (show_segmentation_boundaries + 1) %2;
      break;
      
    case 3: // ctrl-C
      color_scheme = ((color_scheme + 1) % NUM_COLOR_SCHEMES);
      if (color_scheme == PICK_COLOR_SCHEME) color_scheme = RGB_COLOR_SCHEME;
      break;
      
    case 4:   // ctrl-D
    case 8:   // BACKSPACE
    case 127: // DELETE
      RemoveAssignment();
      break;
      
    case 9: // ctrl-I
      show_image = ((show_image + 1) % 2);
      break;
      
    case 12: // ctrl-L
      show_labels = ((show_labels + 1) % 2);
      break;

    case 13: // ctrl-M
      MergeSelectedSegmentations();
      break;

    case 17: // ctrl-Q
      GLUTStop();
      break;

    case 18: // ctrl-R
      show_channels_in_rgb = (show_channels_in_rgb + 1) % 2;
      break;

    case 19: // ctrl-S
      if (!output_keystroke_name) output_keystroke_name = "keystrokes.txt";
      WriteKeystrokes(scene, output_keystroke_name);
      break;
      
    case 26: // ctrl-Z
      // Remove last keystroke
      if (!keystrokes.IsEmpty()) {
        Keystroke *keystroke = keystrokes.Tail();
        keystrokes.RemoveTail();
        delete keystroke;
      }
      break;

    case 32: // Space bar
      display_type = (display_type + 1) % 2;
      viewing_angle = -0.25 * RN_PI;
      break;
      
    case 27: // ESCAPE
      EmptySelections();
      break;

    case '(':
    case ')':
    case '_':
    case '-':    
    case '+':
    case '=':    
      if (selected_channel) {
        RNScalar center = viewing_range.Mid();
        RNScalar radius = viewing_range.Radius();
        if (key == '(') radius *= 0.95;
        else if (key == ')') radius *= 1.05;
        else if (key == '_') radius *= 0.95;
        else if (key == '+') radius *= 1.05;
        else if (key == '-') center -= 0.1 * viewing_range.Radius();
        else if (key == '=') center += 0.1 * viewing_range.Radius();
        if (radius < 0.001) radius = 0.001;
        RNScalar min_value = center - radius;
        RNScalar max_value = center + radius;
        viewing_range.Reset(min_value, max_value);
        glutPostRedisplay();
      }
      break;
    }
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
  if ( glutGetWindow() != GLUTwindow ) 
    glutSetWindow(GLUTwindow);  

  // Redraw
  glutPostRedisplay();
}

#endif



void GLUTInit(int *argc, char **argv)
{
  // Set window dimensions
  GLUTwindow_width = 1200;
  GLUTwindow_height = 600;
  if (selected_streetside) {
    if (selected_channel) {
      // RNScalar aspect = (RNScalar) selected_channel->height / (RNScalar) selected_channel->width;
      // GLUTwindow_height = aspect * GLUTwindow_width;
    }
  }

  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 10);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
  GLUTwindow = glutCreateWindow("pfmview");
  glutSetWindowTitle("gsvlabel");

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
  atexit(AtExit);

  // Remember that glut is initialized
  GLUTinitialized = 1;
}



void GLUTMainLoop(void)
{
  // Run main loop -- never returns 
  glutMainLoop();
}


 
////////////////////////////////////////////////////////////////////////
// Argument parsing
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-interactive")) { interactive = 1; }
      else if (!strcmp(*argv, "-input_labels")) { argc--; argv++; input_labels_name = *argv; }
      else if (!strcmp(*argv, "-input_segmentation")) { argc--; argv++; input_segmentation_name = *argv; }
      else if (!strcmp(*argv, "-input_ground_truth")) { argc--; argv++; input_ground_truth_name = *argv; }
      else if (!strcmp(*argv, "-input_keystrokes")) { argc--; argv++; input_keystroke_name = *argv; }
      else if (!strcmp(*argv, "-output_keystrokes")) { argc--; argv++; output_keystroke_name = *argv; }
      else if (!strcmp(*argv, "-output_scan_rectangles")) { argc--; argv++; output_scan_rectangle_name = *argv; }
      else if (!strcmp(*argv, "-output_image_rectangles")) { argc--; argv++; output_image_rectangle_name = *argv; }
      else if (!strcmp(*argv, "-output_keystrokes")) { argc--; argv++; output_keystroke_name = *argv; }
      else if (!strcmp(*argv, "-parameterization")) { argv++; argc--; parameterization = *argv; }
      else if (!strcmp(*argv, "-DA")) { parameterization = "DA"; }
      else if (!strcmp(*argv, "-SA")) { parameterization = "SA"; }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else channel_names.Insert(*argv);
      argv++; argc--;
    }
  }

  // Check gsv scene filename
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsvlabel gsvfile [options].\n");
    return 0;
  }

  // Check channel names
  if (channel_names.IsEmpty()) channel_names.Insert("ViewpointDepth");

  // Check if should be interactive
  if (!output_scan_rectangle_name && !output_image_rectangle_name) interactive = 1;

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read scene
  scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Read keystrokes
  if (input_keystroke_name) {
    if (!ReadKeystrokes(scene, input_keystroke_name)) exit(-1);
  }

  // Select first streetside
  SelectStreetSide(0);

  // Check if should start interactive interface
  if (interactive) {
    // Initialize GLUT
    GLUTInit(&argc, argv);

    // Run GLUT interface
    GLUTMainLoop();
  }

  // Writes output files
  AtExit();

  // Return success 
  return 0;
}




////////////////////////////////////////////////////////////////////////
// StreetSide functions
////////////////////////////////////////////////////////////////////////

StreetSide::
StreetSide(GSVScan *scan) 
  : scan(scan),
    parse(NULL),
    scanline_grid(NULL),
    point_backpointers(NULL),
    segmentation_backpointers(NULL),
    channels(), 
    index(-1) 
{
}



StreetSide::
~StreetSide(void) 
{
  // Delete everything
  Release();
}




int StreetSide::
Read(void) 
{
  // Get useful variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  GSVScene *scene = run->Scene();
  const char *cache_directory = scene->CacheDataDirectoryName();
  char filename[4096];

  // Read scanline grid
  scanline_grid = new R2Grid();
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_Scanline.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!scanline_grid->Read(filename)) return 0;

  // Read position grids
  R2Grid position_x_grid, position_y_grid, position_z_grid;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_PositionX.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!position_x_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_PositionY.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!position_y_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_PositionZ.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!position_z_grid.Read(filename)) return 0;

  // Read normal grids
  R2Grid normal_x_grid, normal_y_grid, normal_z_grid;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_NormalX.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!normal_x_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_NormalY.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!normal_y_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_NormalZ.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!normal_z_grid.Read(filename)) return 0;

  // Read color grids
  R2Grid color_red_grid, color_green_grid, color_blue_grid;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_ColorRed.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!color_red_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_ColorGreen.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!color_green_grid.Read(filename)) return 0;
  sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_ColorBlue.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!color_blue_grid.Read(filename)) return 0;

  // Read parse grids
  R2Grid parse_detections_grid, parse_labels_grid;
  sprintf(filename, "%s/parse/%s/%02d_%02d_%s_ParseDetections.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (RNFileExists(filename)) {
    sprintf(filename, "%s/parse/%s/%02d_%02d_%s_ParseDetections.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
    if (!parse_detections_grid.Read(filename)) return 0;
    sprintf(filename, "%s/parse/%s/%02d_%02d_%s_ParseLabels.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
    if (!parse_labels_grid.Read(filename)) return 0;
  }
  else {
    sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_%s.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization, input_segmentation_name);
    if (!parse_detections_grid.Read(filename)) return 0;
    parse_labels_grid = parse_detections_grid;  
    parse_labels_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  }

  // Read channel grids
  for (int i = 0; i < channel_names.NEntries(); i++) {
    const char *channel_name = channel_names[i];
    Channel *channel = new Channel(this, parameterization, channel_name);
    if (!channel->Read()) return 0;
    channel->index = channels.NEntries();
    channels.Insert(channel);
  }

  // Initialize backpointers
  point_backpointers = new ObjectPoint * [ scanline_grid->NEntries() ];
  segmentation_backpointers = new ObjectSegmentation * [ scanline_grid->NEntries() ];
  for (int i = 0; i < scanline_grid->NEntries(); i++) point_backpointers[i] = NULL;
  for (int i = 0; i < scanline_grid->NEntries(); i++) segmentation_backpointers[i] = NULL;

  // Create parse
  parse = new ObjectParse();
  sprintf(filename, "%s/parse/%s/%02d_%02d_%s_Parse.ssc", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  printf("%s\n", filename);
  if (RNFileExists(filename)) { if (!parse->ReadFile(filename)) return 0; }
  else if (input_labels_name) { if (!parse->ReadFile(input_labels_name)) return 0; }

  // Create parse detections and segmentations
  assert(parse_detections_grid.NEntries() == scanline_grid->NEntries());
  for (int i = 0; i < parse_detections_grid.NEntries(); i++) {
    // Get position
    RNScalar px = position_x_grid.GridValue(i);
    if (px == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar py = position_y_grid.GridValue(i);
    if (py == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar pz = position_z_grid.GridValue(i);
    if (pz == R2_GRID_UNKNOWN_VALUE) continue;
    R3Point position(px, py, pz);

    // Get normal
    R3Vector normal(0, 0, 0);
    RNScalar nx = normal_x_grid.GridValue(i);
    if (nx != R2_GRID_UNKNOWN_VALUE) {
      RNScalar ny = normal_y_grid.GridValue(i);
      if (ny != R2_GRID_UNKNOWN_VALUE) {
        RNScalar nz = normal_z_grid.GridValue(i);
        if (nz != R2_GRID_UNKNOWN_VALUE) {
          normal.Reset(nx, ny, nz);
        }
      }
    }

    // Get color
    RNRgb color(0.5, 0.5, 0.5);
    RNScalar red = color_red_grid.GridValue(i);
    if (red != R2_GRID_UNKNOWN_VALUE) {
      RNScalar green = color_green_grid.GridValue(i);
      if (green != R2_GRID_UNKNOWN_VALUE) {
        RNScalar blue = color_blue_grid.GridValue(i);
        if (blue != R2_GRID_UNKNOWN_VALUE) {
          color.Reset(red, green, blue);
        }
      }
    }

    // Get detection index
    RNScalar detection_value = parse_detections_grid.GridValue(i);
    if (detection_value == R2_GRID_UNKNOWN_VALUE) continue;
    int detection_index = (int) (detection_value + 0.5);

    // Create detection and segmentation (if necessary)
    while (parse->NDetections() <= detection_index) {
      ObjectDetection *detection = new ObjectDetection();
      parse->InsertDetection(detection);
      ObjectSegmentation *segmentation = new ObjectSegmentation();
      detection->InsertSegmentation(segmentation);
    }

    // Get detection and segmentation
    ObjectDetection *detection = parse->Detection(detection_index);
    assert(detection->NSegmentations() > 0);
    ObjectSegmentation *segmentation = detection->Segmentation(0);

    // Create point
    ObjectPoint local_point(position, normal, color); local_point.SetValue(i);
    segmentation->pointset.InsertPoint(local_point); // THIS IS A HACK
    ObjectPoint *point = segmentation->pointset.points.Tail(); // THIS IS A HACK
    point_backpointers[i] = point;
    segmentation_backpointers[i] = segmentation;

    // Assign descriptor
    if (channels.NEntries() > 0) {
      // Compute descriptor
      ObjectDescriptor& descriptor = point->descriptor;
      descriptor.Reset(NULL, channels.NEntries());
      for (int j = 0; j < channels.NEntries(); j++) {
        Channel *channel = channels.Kth(j);
        for (int k = 0; k < 1; k++) {
          R2Grid *grid = channel->grids.Kth(k);
          RNScalar value = grid->GridValue(i);
          if (value == R2_GRID_UNKNOWN_VALUE) value = 0;
          RNScalar scale = (channel->stddev > 0) ? 1.0 / channel->stddev : 1.0;
          descriptor.SetValue(j, scale * (value - channel->mean));
        }
      }
    }

    // Assign label if there is one
    RNScalar label_value = parse_labels_grid.GridValue(i);
    if (label_value != R2_GRID_UNKNOWN_VALUE) {
      if (segmentation->NAssignments() == 0) {
        // Get label index
        int label_index = (int) (label_value + 0.5);

        // Create label and model (if necessary)
        while (parse->NLabels() <= label_index) {
          ObjectLabel *label = new ObjectLabel();
          parse->InsertLabel(label);
          ObjectModel *model = new ObjectModel();
          label->InsertModel(model);
        }

        // Get label and model
        ObjectLabel *label = parse->Label(label_index);
        assert(label->NModels() > 0);
        ObjectModel *model = label->Model(0);

        // Create assignment
        ObjectAssignment *assignment = new ObjectAssignment(segmentation, model, R3identity_affine, 1.0, TRUE);
        segmentation->InsertAssignment(assignment);
        model->InsertAssignment(assignment);
      }
    }
  }

  // Select channel
  assert(selected_channel_index >= 0);
  assert(selected_channel_index < channels.NEntries());
  SelectChannel(selected_channel_index);

  // Return success
  return 1;
}



int StreetSide::
Write(void)
{
  // Get convenient variables
  if (!parse) return 0;
  if (!scanline_grid) return 0;
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  GSVScene *scene = run->Scene();
  const char *cache_directory = scene->CacheDataDirectoryName();
  char filename[4096];

  // Initialize grids
  R2Grid parse_detections_grid(*scanline_grid);
  parse_detections_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid parse_labels_grid(parse_detections_grid);

  // Fill grids
  for (int detection_index = 0; detection_index < parse->NDetections(); detection_index++) {
    ObjectDetection *detection = parse->Detection(detection_index);
    if (detection->NSegmentations() == 0) continue;
    ObjectSegmentation *segmentation = detection->Segmentation(0);
    const ObjectPointSet& pointset = segmentation->PointSet();
    ObjectLabel *label = segmentation->GroundTruthLabel();
    RNScalar label_index = (label) ? label->ParseIndex() : R2_GRID_UNKNOWN_VALUE;
    for (int k = 0; k < pointset.NPoints(); k++) {
      ObjectPoint *point = pointset.Point(k);
      int grid_index = (int) (point->Value() + 0.5);
      parse_detections_grid.SetGridValue(grid_index, detection_index);
      parse_labels_grid.SetGridValue(grid_index, label_index);
    }
  }

  // Create directory
  sprintf(filename, "mkdir -p %s/parse/%s", cache_directory, run->Name());
  system(filename);

  // Write grids
  sprintf(filename, "%s/parse/%s/%02d_%02d_%s_ParseDetections.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!parse_detections_grid.Write(filename)) return 0;
  sprintf(filename, "%s/parse/%s/%02d_%02d_%s_ParseLabels.grd", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!parse_labels_grid.Write(filename)) return 0;

  // Write parse
  sprintf(filename, "%s/parse/%s/%02d_%02d_%s_Parse.ssc", cache_directory, run->Name(), segment_index, scan_index, parameterization);
  if (!parse->WriteFile(filename)) return 0;

  // Return success
  return 1;
}



void StreetSide::
Release(void) 
{
  // Write parse
  Write();

  // Delete parse
  if (parse) { delete parse; parse = NULL; }

  // Delete grids
  if (scanline_grid) { delete scanline_grid; scanline_grid = NULL; }
  if (point_backpointers) { delete [] point_backpointers; point_backpointers = NULL; }
  if (segmentation_backpointers) { delete [] segmentation_backpointers; segmentation_backpointers = NULL; }

  // Delete channels
  for (int i = 0; i < channels.NEntries(); i++) delete channels[i];
  channels.Empty();
}



void StreetSide::
Draw(int color_scheme) const
{
#if 1
  // Draw selected channel
  if (display_type == 0) DrawGrid(color_scheme);
  else if (display_type == 1) DrawPoints(color_scheme);
#else
  int saved_window_height = GLUTwindow_height;
  glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  DrawGrid(color_scheme);
  glViewport(0, GLUTwindow_height, GLUTwindow_width, GLUTwindow_height);
  DrawPoints(color_scheme);
#endif
  if (show_image) DrawImage(color_scheme);
}



void StreetSide::
DrawImage(int color_scheme) const
{
  // Get channel
  if (!scan) return;

  // Get scanline
  int xcenter = viewing_window.XCenter();
  if (xcenter < 1) xcenter = 1;
  if (xcenter > scanline_grid->XResolution()-1) xcenter = scanline_grid->XResolution()-1;
  int scanline_index = scanline_grid->GridValue(xcenter, 0);
  if (scanline_index == R2_GRID_UNKNOWN_VALUE) return;
  GSVScanline *scanline = scan->Scanline(scanline_index);
  if (!scanline) return;

  // Get image
  RNScalar timestamp = scanline->Timestamp();
  GSVSegment *segment = scan->Segment();
  int tapestry_index = (scan->SegmentIndex() == 0) ? 6 : 2;
  GSVTapestry *tapestry = segment->Tapestry(tapestry_index);
  GSVImage *image = tapestry->FindImageBeforeTimestamp(timestamp);

  // Draw image
  DrawGSVImage(image);
}



void StreetSide::
DrawPoints(int color_scheme) const
{
  // Check position grids
  if (!scan) return;

  // Get window (purposely 4X too big)
  int width = scanline_grid->XResolution();
  int height = scanline_grid->YResolution();
  int xmin = (int) (viewing_window.XCenter() - 2*viewing_window.XLength());
  int xmax = (int) (viewing_window.XCenter() + 2*viewing_window.XLength());
  if (xmin < 0) xmin = 0;
  if (xmax < 0) xmax = 0;
  if (xmin > width-1) xmin = width-1;
  if (xmax > width-1) xmax = width-1;

  // Set projection and modelview matrices
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  PointViewer().Camera().Load();

  // Set pick tolerance
  int pick_tolerance = 10; // (in pixels)
  int point_size = (int) (1.0 * GLUTwindow_width / viewing_window.XLength());
  if (color_scheme == PICK_COLOR_SCHEME) point_size += pick_tolerance;

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);

  // Draw selections
  if (color_scheme != PICK_COLOR_SCHEME) {
    glDisable(GL_DEPTH_TEST);
    glPointSize(1.5 * point_size);
    RNLoadRgb(RNyellow_rgb);
    glBegin(GL_POINTS);
    for (int i = 0; i < selected_segmentations.NEntries(); i++) {
      ObjectSegmentation *segmentation = selected_segmentations.Kth(i);
      const ObjectPointSet& pointset = segmentation->PointSet();
      for (int j = 0; j < pointset.NPoints(); j++) {
        ObjectPoint *point = pointset.Point(j);
        R3LoadPoint(point->Position());
      }
    }
    glEnd();
    glPointSize(1);
    glEnable(GL_DEPTH_TEST);
  }

  // Draw points
  glPointSize(point_size);
  glBegin(GL_POINTS); 
  for (int j = 0; j < height; j++) {
    for (int i = xmin; i <= xmax; i++) {
      int grid_index;
      scanline_grid->IndicesToIndex(i, j, grid_index);
      ObjectPoint *point = point_backpointers[grid_index];
      if (!point) continue;
      if (LoadColor(this, selected_channel, i, j, color_scheme)) {
        R3LoadPoint(point->Position());
      }
    }
  }
  glEnd();
  glPointSize(1);

  // Reset projection and modelview matrices
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Disable depth test
  glDisable(GL_DEPTH_TEST);
}



void StreetSide::
DrawGrid(int color_scheme) const
{
  // Check channel
  if (!selected_channel) return;

  // Get window
  int width = scanline_grid->XResolution();
  int height = scanline_grid->YResolution();
  int xmin = (viewing_window.XMin() > 1) ? viewing_window.XMin() : 1;
  int ymin = (viewing_window.YMin() > 1) ? viewing_window.YMin() : 1;
  int xmax = (viewing_window.XMax()+1 < width-1) ? viewing_window.XMax()+1 : width-1;
  int ymax = (viewing_window.YMax()+1 < height-1) ? viewing_window.YMax()+1 : height-1;

  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(viewing_window.XMin(), viewing_window.XMax(), viewing_window.YMin(), viewing_window.YMax()); 

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Set point size
  int point_size = (int) (1.75 * GLUTwindow_width / viewing_window.XLength());

  // Draw grid
  glPointSize(point_size);
  glBegin(GL_POINTS);
  for (int j = ymin; j <= ymax; j++) {
    for (int i = xmin; i <= xmax; i++) {
      if (LoadColor(this, selected_channel, i, j, color_scheme)) {
        glVertex2i(i, j);
      }
    }
  }
  glEnd();

  // Draw other stuff
  if (color_scheme != PICK_COLOR_SCHEME) {
    // Draw segmentation boundaries
    if (show_segmentation_boundaries) {
      glPointSize(0.5*point_size);
      RNLoadRgb(1.0, 0.0, 0.0);
      glBegin(GL_POINTS);
      for (int j = ymin; j <= ymax; j++) {
        for (int i = xmin; i <= xmax; i++) {
          // Check if segmentation boundary
          int grid_index;
          scanline_grid->IndicesToIndex(i, j, grid_index);
          ObjectSegmentation *segmentation = segmentation_backpointers[grid_index];
          for (int ix = i+1; ix <= i+1; ix+=2) {
            if (ix < 0) continue;
            if (ix >= scanline_grid->XResolution()) continue;
            scanline_grid->IndicesToIndex(ix, j, grid_index);
            if (segmentation_backpointers[grid_index] != segmentation) {
              glVertex2i(i, j);
              break;
            }
          }
          for (int iy = j+1; iy <= j+1; iy+=2) {
            if (iy < 0) continue;
            if (iy >= scanline_grid->YResolution()) continue;
            scanline_grid->IndicesToIndex(i, iy, grid_index);
            if (segmentation_backpointers[grid_index] != segmentation) {
              glVertex2i(i, j);
              break;
            }
          }
        }
      }
      glEnd();
    }

    // Draw selections
    glPointSize(1.5*point_size);
    RNLoadRgb(RNyellow_rgb);
    glBegin(GL_POINTS);
    for (int i = 0; i < selected_segmentations.NEntries(); i++) {
      ObjectSegmentation *segmentation = selected_segmentations.Kth(i);
      const ObjectPointSet& pointset = segmentation->PointSet();
      for (int j = 0; j < pointset.NPoints(); j++) {
        int ix, iy;
        ObjectPoint *point = pointset.Point(j);
        int grid_index = (int) (point->Value() + 0.5);
        scanline_grid->IndexToIndices(grid_index, ix, iy);
        R3LoadPoint(ix, iy, 0.0);
      }
    }
    glEnd();
  }

  // Restore point size
  glPointSize(1);

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}



R3Viewer StreetSide::
PointViewer(void) const
{
  // Used to remember last viewer
  static R3Viewer last_viewer(R3default_camera, R2default_viewport);

  // Return viewer for points display type
  R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  int width = scanline_grid->XResolution();
  int xcenter = viewing_window.XCenter();
  if (xcenter < 1) xcenter = 1;
  if (xcenter > width-1) xcenter = width-1;
  int scanline_index = scanline_grid->GridValue(xcenter, 0);
  if (scanline_index == R2_GRID_UNKNOWN_VALUE) return last_viewer;
  GSVScanline *scanline = scan->Scanline(scanline_index);
  const GSVPose& pose = scanline->Pose();

  // Get rotated viewing vectors
  R3Vector towards = pose.Towards();
  towards.Rotate(pose.Right(), viewing_angle);
  towards.Normalize();
  R3Vector up = pose.Right() % towards;
  up.Normalize();

  // Get lookat point
  RNLength d = (viewing_window.YLength()) / 4;
  int ycenter = viewing_window.YCenter();
  R3Point lookat = pose.Viewpoint() + 5.0*towards + 0.1*ycenter*up;
  R3Point viewpoint = lookat - d*towards;

  // Return viewer
  R3Camera camera(viewpoint, towards, up, 0.5, 0.5, 0.1, 10000.0);
  R3Viewer viewer(camera, viewport);
  last_viewer = viewer;
  return viewer;
}



////////////////////////////////////////////////////////////////////////
// Channel functions
////////////////////////////////////////////////////////////////////////

Channel::
Channel(StreetSide *streetside,  const char *parameterization, const char *name, const char *extension)
  : streetside(streetside),
    parameterization(parameterization ? strdup(parameterization) : NULL),
    name(name ? strdup(name) : NULL),
    extension(extension ? strdup(extension) : NULL),
    grids(),
    mean(0),
    stddev(0),
    minimum(0),
    maximum(0),
    width(0),
    height(0),
    index(-1)
{
}


Channel::
~Channel(void)
{
  // Release everything
  Release();
}



int Channel::
Read(void)
{
  // Check parameters
  if (!parameterization) return 0;
  if (!name) return 0;
  if (!extension) return 0;

  // Get filename
  char filename[4096];
  GSVScan *scan = streetside->scan;
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  GSVScene *scene = run->Scene();
  const char *cache_directory = scene->CacheDataDirectoryName();

  // Check extension
  if (!strcmp(extension, "grd") || !strcmp(extension, "pfm")) {
    // Determine number of components by name
    RNArray<const char *> suffix;       
    if (!strcmp(name, "Color")) { suffix.Insert("Red"); suffix.Insert("Green"); suffix.Insert("Blue"); }
    else if (!strcmp(name, "Position")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Normal")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Normal")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Viewpoint")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Towards")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Up")) { suffix.Insert("X"); suffix.Insert("Y"); suffix.Insert("Z"); }
    else if (!strcmp(name, "Lambda")) { suffix.Insert("1"); suffix.Insert("2"); suffix.Insert("3"); }
    else if (!strcmp(name, "LargePlane")) { suffix.Insert("A"); suffix.Insert("B"); suffix.Insert("C"); }
    else if (!strcmp(name, "SmallPlane")) { suffix.Insert("A"); suffix.Insert("B"); suffix.Insert("C"); }
    else { suffix.Insert(""); }
    for (int i = 0; i < suffix.NEntries(); i++) {
      char filename[4096];
      R2Grid *grid = new R2Grid();
      sprintf(filename, "%s/laser_images/%s/%02d_%02d_%s_%s%s.%s", cache_directory, run->Name(), segment_index, scan_index, parameterization, name, suffix[i], extension);
      if (!grid->Read(filename)) return 0;
      grids.Insert(grid);
      width = grid->XResolution();
      height = grid->YResolution();
    }
  }
  else {
    // Read image from file
    R2Image *image = new R2Image();
    if (!image->Read(filename)) return 0;
    width = image->Width();
    height = image->Height();
    for (int c = 0; c < 3; c++) {
      R2Grid *grid = new R2Grid(width, height);
      grids.Insert(grid);
      for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
          RNRgb color = image->PixelRGB(i, j);
          grid->SetGridValue(i, j, color[c]);
        }
      }
    }
  }

  // Update statistics
  if (grids.NEntries() > 0) {
    R2Grid tmp(*(grids[0]));
    int cardinality = tmp.Cardinality();
    RNInterval range = tmp.Range();
    minimum = range.Min();
    maximum = range.Max();
    mean = tmp.Mean();
    tmp.Subtract(mean);
    RNScalar l2norm = tmp.L2Norm();
    stddev = (cardinality > 0) ? sqrt(l2norm * l2norm / cardinality) : 0.0;
    printf("%g %g\n", mean, stddev);
  }

  // Return success
  return 1;
}



void Channel::
Release(void)
{
  // Delete everything
  if (parameterization) {free(parameterization); parameterization = NULL; }
  if (name) { free(name); name = NULL; }
  if (extension) { free(extension); extension = NULL; }
  for (int i = 0; i < grids.NEntries(); i++) delete grids[i]; 
  grids.Empty();
}



////////////////////////////////////////////////////////////////////////
// Classifier functions
////////////////////////////////////////////////////////////////////////

Classifier::
Classifier(void) 
{
}



Classifier::
~Classifier(void) 
{
}



int Classifier::
Read(const char *filename)
{
  return 0;
}



int Classifier::
Write(const char *filename) const
{
  return 0;
}



int Classifier::
Train(const RNArray<StreetSide *>& streetsides, const RNArray<const char *>& classifier_feature_names) 
{
  return 0;
}



int Classifier::
Predict(const RNArray<StreetSide *>& streetsides, const RNArray<const char *>& classifier_feature_names, const RNArray<const char *>& crf_feature_names) 
{
  return 0;
}







