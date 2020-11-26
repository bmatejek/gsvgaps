/////////////////////////////////////////////////////////////////////////
// Include file for object label class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class ObjectLabel {
public:
  ////////////////////////////////////////
  // Constructor/destructor functions
  ////////////////////////////////////////

  // Constructor destructor
  ObjectLabel(const char *name = NULL, const R3Point& origin = R3unknown_point, 
    RNScalar max_radius = RN_UNKNOWN, RNScalar max_height = RN_UNKNOWN,
    int assignment_keystroke = 0, const RNRgb& color = RNRgb(0.5, 0.5, 0.5));
  ~ObjectLabel(void);


  ////////////////////////////////////////
  // Access functions
  ////////////////////////////////////////

  // Parse access
  ObjectParse *Parse(void) const;
  int ParseIndex(void) const;

  // Model access
  int NModels(void) const;
  ObjectModel *Model(int k) const;

  // Property access
  const char *Name(void) const;
  const R3Point& Origin(void) const;
  RNLength MaxRadius(void) const;
  RNLength MaxHeight(void) const;
  int AssignmentKeystroke(void) const;
  const RNRgb& Color(void) const;

  ////////////////////////////////////////
  // Manipulation functions
  ////////////////////////////////////////

  // Model manipulation
  void InsertModel(ObjectModel *model);
  void RemoveModel(ObjectModel *model);

  // Property manipulation
  void SetName(const char *name);
  void SetOrigin(const R3Point& origin);
  void SetMaxRadius(RNLength radius);
  void SetMaxHeight(RNLength height);
  void SetAssignmentKeystroke(int key);
  void SetColor(const RNRgb& color);

  ////////////////////////////////////////
  // Internal stuff
  ////////////////////////////////////////

public:
  friend class ObjectParse;
  ObjectParse *parse;
  int parse_index;
  RNArray<ObjectModel *> models;
  char *name;
  R3Point origin;
  RNLength max_radius;
  RNLength max_height;
  int assignment_keystroke;
  RNRgb color;
};



////////////////////////////////////////////////////////////////////////
// Extern functions (used for sorting)
////////////////////////////////////////////////////////////////////////

extern int CompareObjectLabels(const void *data1, const void *data2);



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline ObjectParse *ObjectLabel::
Parse(void) const
{
  // Return parse
  return parse;
}



inline int ObjectLabel::
ParseIndex(void) const
{
  // Return parse index
  return parse_index;
}



inline int ObjectLabel::
NModels(void) const
{
  // Return number of models
  return models.NEntries();
}



inline ObjectModel *ObjectLabel::
Model(int k) const
{
  // Return kth model
  return models.Kth(k);
}



inline const char *ObjectLabel::
Name(void) const
{
  // Return name
  return name;
}



inline const R3Point& ObjectLabel::
Origin(void) const
{
  // Return origin
  return origin;
}



inline RNLength ObjectLabel::
MaxRadius(void) const
{
  // Return maximum radius
  return max_radius;
}



inline RNLength ObjectLabel::
MaxHeight(void) const
{
  // Return maximum height
  return max_height;
}



inline int ObjectLabel::
AssignmentKeystroke(void) const
{
  // Return assignment keystroke
  return assignment_keystroke;
}


inline const RNRgb& ObjectLabel::
Color(void) const
{
  // Return color
  return color;
}


