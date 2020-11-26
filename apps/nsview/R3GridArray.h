////////////////////////////////////////////////////////////////////////
// R3Grid array class definition
////////////////////////////////////////////////////////////////////////

struct R3GridArray {
  // Constructors
  R3GridArray(void);
  R3GridArray(const R3GridArray& grids);
  ~R3GridArray(void);

  // Access
  int NGrids(void) const;
  R3Grid *Grid(int k) const;

  // Grid manipulation
  void Empty(void);
  R3GridArray& operator=(const R3GridArray& grids);
  void Insert(const R3GridArray& array);
  void Insert(R3Grid *grid);
  void Remove(R3Grid *grid);

  // Read/write functions
  int ReadFile(const char *filename);
  int ReadMetaFile(const char *filename);
  int ReadGridFile(const char *filename);
  int ReadRawFile(const char *filename, const char *data_type, int xres, int yres, int zres, int ngrids);
  int WriteFile(const char *filename) const;
  int WriteMetaFile(const char *filename) const;
  int WriteGridFile(const char *filename) const;
  int WriteRawFile(const char *filename, const char *data_type = "Float32") const;

public:
  RNArray<struct R3Grid *> grids;
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline int R3GridArray::
NGrids(void) const
{
  return grids.NEntries();
}



inline R3Grid *R3GridArray::
Grid(int k) const
{
  return grids.Kth(k);
}



