#ifndef GENERATE_GRID_H
#define GENERATE_GRID_H

////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

// class declaration to suppress errors
class ObjectOfInterest;

#include <iostream>
#include <fstream>
#include "helper.h"
#include "VOCBox.h"

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

class ObjectOfInterest {
 public:
  // constructor
  ObjectOfInterest(char *object_name, int index);

  // accessors/manipulators
  void SetHeight(double height);

  const char *ObjectName(void) const;
  const double Height(void) const;
  const int Index(void) const;

 private:
  // instance variables
  char *object_name;
  double height;
  int index;
};

class GenerateGrid {
 public:
  // constructors
  GenerateGrid(GSVScene *scene, GSVPoseOptimization *optimization, bool *show_runs, std::vector<ObjectOfInterest *> objects);
  
  // accessors/manipulators
  const GSVScene *Scene(void) const;
  const std::vector<R3Grid *> Grid(void) const;
  const GSVPoseOptimization *Optimization(void) const;
  const bool *ShowRuns(void) const;
  const std::vector<ObjectOfInterest *> ObjectsOfInterest(void) const;

 private:
  // instance variables
  GSVScene *scene;
  std::vector<R3Grid *> grid;
  GSVPoseOptimization *optimization;
  bool *show_runs;
  std::vector<ObjectOfInterest *> objects_of_interest;
};

#endif // GENERATE_GRID_H
