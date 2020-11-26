////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <fstream>
#include <RNBasics\RNBasics.h>
#include <R2Shapes\R2Shapes.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// Input/output
static char *imname = NULL;
static int print_verbose = 0;

////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage(void) {
	fprintf(stderr, "Usage: txt2grd imname [-v]\n");
	return 0;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
			else if (!strcmp(*argv, "-usage")) { return Usage(); }
			else {

				fprintf(stderr, "Invalid program argument: %s\n", *argv);
				return Usage();
			}
			argv++; argc--;
		}
		else {
			if (!imname) imname = *argv;
			else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
	}

	// Check scene name
	if (!imname) return Usage();

	// Return OK status 
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	RNTime start_time;
	start_time.Read();

	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

	char input_filename[4096];

	ifstream input_class;
	sprintf(input_filename, "output/%s_classID.txt", imname);
	input_class.open(input_filename);
	if (!input_class) { fprintf(stderr, "Failed to read %s\n", input_filename); exit(-1); }

	ifstream input_score;
	sprintf(input_filename, "output/%s_score.txt", imname);
	input_score.open(input_filename);
	if (!input_score) { fprintf(stderr, "Failed to read %s\n", input_filename); exit(-1); }

	ifstream input_object;
	sprintf(input_filename, "output/%s_objectID.txt", imname);
	input_object.open(input_filename);
	if (!input_object) { fprintf(stderr, "Failed to read %s\n", input_filename); exit(-1); }

	string classline;
	string scoreline;
	string objectline;

	int height = 0;
	int width = 0;
	while (getline(input_class, classline) && getline(input_score, scoreline) && getline(input_object, objectline)) {
		width = splitString(classline, ',').size();
		height++;
	}
	input_class.clear();
	input_class.seekg(0);
	input_score.clear();
	input_score.seekg(0);
	input_object.clear();
	input_object.seekg(0);

	R2Grid classgrd = R2Grid(width, height);
	R2Grid scoregrd = R2Grid(width, height);
	R2Grid objectgrd = R2Grid(width, height);

	int y = 0;
	while (getline(input_class, classline) && getline(input_score, scoreline) && getline(input_object, objectline)) {
		vector<string> classes = splitString(classline, ',');
		vector<string> scores = splitString(scoreline, ',');
		vector<string> objects = splitString(objectline, ',');
		
		for (unsigned int x = 0; x < classes.size(); x++) {
			classgrd.SetGridValue(x, height - 1 - y, atoi(classes[x].c_str()));
			scoregrd.SetGridValue(x, height -1 - y, atof(scores[x].c_str()));
			objectgrd.SetGridValue(x, height - 1 - y, atoi(objects[x].c_str()));
		}

		y++;
	}

	char output_filename[4096];
	sprintf(output_filename, "grds/%s_ClassID.grd", imname);
	classgrd.Write(output_filename);
	sprintf(output_filename, "grds/%s_Score.grd", imname);
	scoregrd.Write(output_filename);
	sprintf(output_filename, "grds/%s_ObjectID.grd", imname);
	objectgrd.Write(output_filename);

	if (print_verbose) {
		printf("Created output grids ... \n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
	}

	// Return success 
	return 0;
}


