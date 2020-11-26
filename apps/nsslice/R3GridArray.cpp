////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"
#include "R3GridArray.h"


////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

R3GridArray::
R3GridArray(void)
  : grids()
{
}



R3GridArray::
R3GridArray(const R3GridArray& array)
  : grids()
{
  // Copy grids
  for (int i = 0; i < array.NGrids(); i++) {
    R3Grid *grid2 = array.Grid(i);
    R3Grid *grid = new R3Grid(*grid2);
    Insert(grid);
  }
}



R3GridArray::
~R3GridArray(void)
{
  Empty();
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void R3GridArray::
Empty(void)
{
  // Delete the grids
  for (int i = 0; i < grids.NEntries(); i++) {
    delete grids[i];
  }

  // Empty the array
  grids.Empty();
}



void R3GridArray::
Insert(const R3GridArray& array)
{
  // Copy grids
  for (int i = 0; i < array.NGrids(); i++) {
    R3Grid *grid2 = array.Grid(i);
    R3Grid *grid = new R3Grid(*grid2);
    Insert(grid);
  }
}



void R3GridArray::
Insert(R3Grid *grid)
{
  // Insert grid
  grids.Insert(grid);
}



void R3GridArray::
Remove(R3Grid *grid)
{
  // Remove grid
  grids.Remove(grid);
}



R3GridArray& R3GridArray::
operator=(const R3GridArray& array)
{
  // Empty existing grids
  Empty();

  // Copy grids
  for (int i = 0; i < array.NGrids(); i++) {
    R3Grid *grid2 = array.Grid(i);
    R3Grid *grid = new R3Grid(*grid2);
    Insert(grid);
  }

  // Return this
  return *this;
}



////////////////////////////////////////////////////////////////////////
// Input/output
////////////////////////////////////////////////////////////////////////

int R3GridArray::
ReadFile(const char *filename)
{
  // Check filename
  if (!filename) return 0;

  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .meta)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".meta", 5)) {
    if (!ReadMetaFile(filename)) return 0;
  }
  else {
    if (!ReadGridFile(filename)) return 0;
  }

  // Return success
  return 1;
}



int R3GridArray::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .fet)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".meta", 5)) {
    if (!WriteMetaFile(filename)) return 0;
  }
  else {
    if (!WriteGridFile(filename)) return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Meta file Input/output
////////////////////////////////////////////////////////////////////////

int R3GridArray::
ReadMetaFile(const char *meta_filename) 
{
  // Open meta file
  FILE *fp = fopen(meta_filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open meta file %s\n", meta_filename);
    return 0;
  }

  // Read meta file
  char buffer[4096];
  if (!fgets(buffer, 4096, fp)) {
    fprintf(stderr, "Unable to read meta file %s\n", meta_filename);
    return 0;
  }

  // Parse data type
  char *data_type = strtok(buffer, "(");
  if (!data_type) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return 0; }

  // Parse dimensions
  char *bufferp = strtok(NULL, ",");
  if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return 0; }
  int xres = atoi(bufferp);
  bufferp = strtok(NULL, ",");
  if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return 0; }
  int yres = atoi(bufferp);
  bufferp = strtok(NULL, ",)");
  if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return 0; }
  int zres = atoi(bufferp);
  bufferp = strtok(NULL, ")");
  int ngrids = (bufferp) ? atoi(bufferp) : 1;

  // Close meta file
  fclose(fp);

  // Read the raw file
  char raw_filename[4096];
  strncpy(raw_filename, meta_filename, 4096);
  char *endp = strrchr(raw_filename, '.');
  if (endp) *endp = '\0';
  strncat(raw_filename, ".raw", 4096);
  if (!ReadRawFile(raw_filename, data_type, xres, yres, zres, ngrids)) return 0;

  // Return success
  return 1;
}



int R3GridArray::
WriteMetaFile(const char *meta_filename) const
{
  // Get meta info
  int ngrids = NGrids();
  if (ngrids == 0) return 0;
  int xres = Grid(0)->XResolution();
  if (xres == 0) return 0;
  int yres = Grid(0)->YResolution();
  if (yres == 0) return 0;
  int zres = Grid(0)->ZResolution();
  if (zres == 0) return 0;
  const char *data_type = "Float32";

  // Open meta file
  FILE *fp = fopen(meta_filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open meta file %s\n", meta_filename);
    return 0;
  }

  // Write meta file
  if (NGrids() > 1) fprintf(fp, "%s(%d,%d,%d,%d)\n", data_type, xres, yres, zres, ngrids);
  else fprintf(fp, "%s(%d,%d,%d)\n", data_type, xres, yres, zres);

  // Close meta file
  fclose(fp);

  // Write the raw file
  char raw_filename[4096];
  strncpy(raw_filename, meta_filename, 4096);
  char *endp = strrchr(raw_filename, '.');
  if (endp) *endp = '\0';
  strncat(raw_filename, ".raw", 4096);
  if (!WriteRawFile(raw_filename, data_type)) return 0;

  // Return success
  return 1;
}




////////////////////////////////////////////////////////////////////////
// Raw file Input/output
////////////////////////////////////////////////////////////////////////

static int
ReadRawValue(FILE *fp, int format, RNScalar &a) 
{
  switch (format) {
  case 1: {
    unsigned char value;
    fread(&value, sizeof(unsigned char), 1, fp);
    a = (RNScalar) value;
    break; }

  case 2: {
    unsigned short value;
    fread(&value, sizeof(unsigned short), 1, fp);
    a = (RNScalar) value;
    break; }

  case 4: {
    float value;
    fread(&value, sizeof(float), 1, fp);
    a  = (RNScalar) value;
    break; }

   case 8: {
    double value;
    fread(&value, sizeof(double), 1, fp);
    a = (RNScalar) value;
    break; }

  default:
    fprintf(stderr, "Unrecognized format %d\n", format);
    return 0;
  }

  // Return success
  return 1;
}



static int
WriteRawValue(FILE *fp, int format, RNScalar a) 
{
  switch (format) {
  case 1: {
    unsigned char value = (unsigned char) a;
    fwrite(&value, sizeof(unsigned char), 1, fp);
    break; }

  case 2: {
    unsigned short value = (unsigned short) a;;
    fwrite(&value, sizeof(unsigned short), 1, fp);
    break; }

  case 4: {
    float value = (float) a;
    fwrite(&value, sizeof(float), 1, fp);
    break; }

   case 8: {
    double value = (double) a;
    fwrite(&value, sizeof(double), 1, fp);
    break; }

  default:
    fprintf(stderr, "Unrecognized format: %d\n", format);
    return 0;
  }

  // Return success
  return 1;
}



int R3GridArray::
ReadRawFile(const char *filename, const char *data_type, int xres, int yres, int zres, int ngrids) 
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw file %s\n", filename);
    return 0;
  }

  // Determine data format
  int format = 0;
  if (!strcmp(data_type, "Uint16")) format = 2; 
  if (!strcmp(data_type, "Float32")) format = 4; 

  // Read file
  for (int i = 0; i < ngrids; i++) {
    R3Grid *grid = new R3Grid(xres, yres, zres);
    for (int iz = 0; iz < zres; iz++) {
      for (int iy = 0; iy < yres; iy++) {
        for (int ix = 0; ix < xres; ix++) {
          RNScalar value = 0;
          if (!ReadRawValue(fp, format, value)) return 0;
          grid->SetGridValue(ix, iy, iz, value);
        }
      }
    }
    Insert(grid);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3GridArray::
WriteRawFile(const char *filename, const char *data_type) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw file %s\n", filename);
    return 0;
  }

  // Determine data format
  int format = 0;
  if (!strcmp(data_type, "Uint16")) format = 2; 
  if (!strcmp(data_type, "Float32")) format = 4; 

  // Write file
  for (int i = 0; i < NGrids(); i++) {
    R3Grid *grid = Grid(i);
    for (int iz = 0; iz < grid->ZResolution(); iz++) {
      for (int iy = 0; iy < grid->YResolution(); iy++) {
        for (int ix = 0; ix < grid->XResolution(); ix++) {
          RNScalar value = grid->GridValue(ix, iy, iz);
          if (!WriteRawValue(fp, format, value)) return 0;
        }
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Grid file Input/output
////////////////////////////////////////////////////////////////////////

int R3GridArray::
ReadGridFile(const char *filename) 
{
  // Read one grid 
  R3Grid *grid = new R3Grid();
  if (!grid->ReadFile(filename)) return 0;
  Insert(grid);
  return 1;
}



int R3GridArray::
WriteGridFile(const char *filename) const
{
  // Write first grid 
  if (grids.IsEmpty()) return 0;
  if (!grids[0]->WriteFile(filename)) return 0;
  return 1;
}




