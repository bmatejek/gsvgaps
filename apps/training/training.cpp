////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

// GSV include files
#include "GSV/GSV.h"

// std library include files
#include <vector>
#include <fstream>
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

class ObjectLabel {
public:
    ObjectLabel(string object, string color_hex, RNRgb color, int id);
    
public:
    string object;
    string color_hex;
    RNRgb color;
    int id;
    int occurences;
};

ObjectLabel::ObjectLabel(string object, string color_hex, RNRgb color, int id) :
object(object),
color_hex(color_hex),
color(color),
id(id)
{
    occurences = 0;
}

static int DA_images = 0;
static int SA_images = 0;
static int segment_index = 0;
static int trainval = 0;
static int scan_index = 0;
static int print_verbose = 0;
static vector<char *> run_names = vector<char *>();
static vector<GSVRun *> runs = vector<GSVRun *>();
static vector<ObjectLabel *> objects = vector<ObjectLabel *>();

////////////////////////////////////////////////////////////////////////
// Basic Input/Output
////////////////////////////////////////////////////////////////////////

static GSVScene *ReadScene(const char *filename)
{
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

////////////////////////////////////////////////////////////////////////
// Training functions
////////////////////////////////////////////////////////////////////////

static void CreateTraining(GSVRun *run, char *parametrization) {
    for (int is = 0; is < run->NSegments(); ++is) {
        GSVSegment *segment = run->Segment(is);
        for (int ia = 0; ia < segment->NScans(); ++ia) {
            if (ia == 1) continue;
            // read in the parse file
            R2Grid parse_labels;
            R2Grid viewpoint_distance;
            R2Image color;

            // read in all of the grid files
            char image_name[4096];
            sprintf(image_name, "gsv_data/parse/%s/%02d_%02d_%s_ParseLabels.grd", run->Name(), segment_index, scan_index, parametrization);

            // make sure the parse file exists
            ifstream ifile(image_name);
            if (!ifile) { continue; }


            parse_labels.Read(image_name);
            sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_ViewpointDistance.grd", run->Name(), segment_index, scan_index, parametrization);
            viewpoint_distance.Read(image_name);
            sprintf(image_name, "gsv_data/laser_images/%s/%02d_%02d_%s_Color.bmp", run->Name(), segment_index, scan_index, parametrization);
            color.Read(image_name);

            // get width and height of this segment
            int width = color.Width();
            int height = color.Height();

            // create all of the images
            int total_images = width / 500;
            for (int t = 0; t < total_images; t++) {
                int start = t * 500;
                int end;

                // get the size of this segment
                if (t == total_images - 1) { end = color.Width(); }
                else { end = (t + 1) * 500; }

                // create new images
                R2Image color_image = R2Image(end - start, height);
                R2Image viewpoint_distance_image = R2Image(end - start, height);
                R2Image parse_labels_image = R2Image(end - start, height);

                double maximum = viewpoint_distance.Maximum();

                for (int i = start; i < end; i++) {
                    for (int j = 0; j < height; j++) {
                        color_image.SetPixelRGB(i - start, j, color.PixelRGB(i, j));

                        // get the pixel for the viewpoint distance
                        RNScalar viewpoint_distance_value = viewpoint_distance.GridValue(i, j);
                        if (viewpoint_distance_value == R2_GRID_UNKNOWN_VALUE) {
                            viewpoint_distance_value = 0;
                        }
                        viewpoint_distance_value = viewpoint_distance_value / maximum;
                        RNRgb pixel = RNRgb(viewpoint_distance_value, viewpoint_distance_value, viewpoint_distance_value);
                        viewpoint_distance_image.SetPixelRGB(i - start, j, pixel);

                        // get the pixel color for the label
                        RNScalar label_value = parse_labels.GridValue(i, j);
                        pixel = RNRgb(0.0, 0.0, 0.0);
                        // get color from id
                        for (unsigned int ob = 0; ob < objects.size(); ++ob) {
                            if (objects[ob]->id == (int)round(label_value)) {
                                pixel = objects[ob]->color;
                                objects[ob]->occurences++;
                                break;
                            }
                        }
                        parse_labels_image.SetPixelRGB(i - start, j, pixel);
                    }
                }

                // create the color image
                sprintf(image_name, "trainval/%s_Color/%s_%02d_%02d_%s_%06d_%06d_Color.png", parametrization, run->Name(), segment_index, scan_index, parametrization, start, end);
                color_image.Write(image_name);

                // create the viewpoint image
                sprintf(image_name, "trainval/%s_ViewpointDistance/%s_%02d_%02d_%s_%06d_%06d_ViewpointDistance.png", parametrization, run->Name(), segment_index, scan_index, parametrization, start, end);
                viewpoint_distance_image.Write(image_name);

                // create the parse labels image
                sprintf(image_name, "labels/%s_%02d_%02d_%s_%06d_%06d_ParseLabels.png", run->Name(), segment_index, scan_index, parametrization, start, end);
                parse_labels_image.Write(image_name);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int Usage() {
	fprintf(stderr, "Usage: training runs -v\n");
	return 0;
}

static int ParseArgs(int argc, char **argv) {
	// Parse arguments
	argc--; argv++;
	while (argc > 0) {
		if ((*argv)[0] == '-') {
			if (!strcmp(*argv, "-v")) { print_verbose = 1; }
            else if (!strcmp(*argv, "-SA")) { SA_images = 1; }
            else if (!strcmp(*argv, "-DA")) { DA_images = 1; }
            else if (!strcmp(*argv, "-scan")) { argv++; argc--; scan_index = atoi(*argv); }
            else if (!strcmp(*argv, "-segment")) { argv++; argc--; segment_index = atoi(*argv); }
            else if (!strcmp(*argv, "-trainval")) { trainval = 1; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return Usage(); }
			argv++; argc--;
		}
		else {
            run_names.push_back(*argv);
            argv++; argc--;
		}
	}

    // Return OK status 
	return 1;
}

static string ConvertSingleHex(int col) {
    if (col == 0) return "00";
    else if (col == 1) return "11";
    else if (col == 2) return "22";
    else if (col == 3) return "33";
    else if (col == 4) return "44";
    else if (col == 5) return "55";
    else if (col == 6) return "66";
    else if (col == 7) return "77";
    else if (col == 8) return "88";
    else if (col == 9) return "99";
    else if (col == 10) return "AA";
    else if (col == 11) return "BB";
    else if (col == 12) return "CC";
    else if (col == 13) return "DD";
    else if (col == 14) return "EE";
    else if (col == 15) return "FF";
    else {
        fprintf(stderr, "Error with input to ConvertSingleHex: %d\n", col);
        exit(-1);
    }
}

static string ConvertToHex(int r, int g) {
    return "0x" + ConvertSingleHex(r) + ConvertSingleHex(g) + "FF";
}

static double ConvertFromHex(string col) {
    if (col == "00") return 0;
    else if (col == "11") return 17;
    else if (col == "22") return 34;
    else if (col == "33") return 51;
    else if (col == "44") return 68;
    else if (col == "55") return 85;
    else if (col == "66") return 102;
    else if (col == "77") return 119;
    else if (col == "88") return 136;
    else if (col == "99") return 153;
    else if (col == "AA") return 170;
    else if (col == "BB") return 187;
    else if (col == "CC") return 204;
    else if (col == "DD") return 221;
    else if (col == "EE") return 238;
    else if (col == "FF") return 255;
    else {
        fprintf(stderr, "Error with input to ConvertFromHex: %s\n", col.c_str());
        exit(-1);
    }
}

static RNRgb PixelColor(string hex) {
    string r = hex.substr(2, 2);
    string g = hex.substr(4, 2);
    string b = hex.substr(6, 2);
    return RNRgb(ConvertFromHex(r) / 255.0, ConvertFromHex(g) / 255.0, ConvertFromHex(b) / 255.0);
}

static void CreateSetFile() {
    for (unsigned int ir = 0; ir < runs.size(); ++ir) {
        GSVRun *run = runs[ir];
        for (int is = 0; is < run->NSegments(); ++is) {
            GSVSegment *segment = run->Segment(is);
            for (int ia = 0; ia < segment->NScans(); ++ia) {
                if (ia == 1) continue;
                char parse_filename[4096];
                sprintf(parse_filename, "gsv_data/parse/%s/%02d_%02d_SA_Parse.ssc", run->Name(), is, ia);
                ifstream ssc_file;
                ssc_file.open(parse_filename);
                if (!ssc_file.is_open()) { fprintf(stderr, "Failed to read %s\n", parse_filename); continue; }
                // read in all of the lines
                string line;
                int id = 0;
                while (getline(ssc_file, line)) {
                    vector<string> parse_file = splitString(line, ' ');
                    if (strcmp("L", parse_file[0].c_str()) && id != 0) break;
                    string name = parse_file[1].c_str();
                    if (name == "PARSE") continue;
                    
                    bool found = false;
                    //printf("%s\n", line.c_str());
                    for (unsigned int ob = 0; ob < objects.size(); ob++) {
                        if (!name.compare(objects[ob]->object) && id != objects[ob]->id) {
                            fprintf(stderr, "Error in reading parse files %d: %s and %d: %s\n", id, name.c_str(), objects[ob]->id, objects[ob]->object.c_str());
                            exit(-1);
                        }
                        else if (!name.compare(objects[ob]->object) && id == objects[ob]->id) {
                            found = true;
                        }
                    }

                    if (!found) {
                        int r = id % 16;
                        int g = 1 + id / 16;
                        if (g > 16) { fprintf(stderr, "Too many objects\n"); }
                        string hex_number = ConvertToHex(r, g);
                        ObjectLabel *obj = new ObjectLabel(string(name).c_str(), hex_number, PixelColor(hex_number), id);
                        objects.push_back(obj);
                    }
                    ++id;
                }
                ssc_file.close();
            }
        }
    }
}

static void SaveSetFile() {
    ofstream dataset_config;
    dataset_config.open("dataset_config.set");


    for (unsigned int ob = 0; ob < objects.size(); ++ob) {
        if (objects[ob]->occurences == 0) {
            objects.erase(objects.begin() + ob);
            ob--;
        }
    }

    for (unsigned int ob = 0; ob < objects.size(); ++ob) {
        printf("%d %s %d\n", objects[ob]->id, objects[ob]->object.c_str(), objects[ob]->occurences);
    }


    // output classes
    dataset_config << "classes=";
    dataset_config << objects.size();
    dataset_config << '\n';
    for (unsigned int ob = 0; ob < objects.size(); ++ob) {
        dataset_config << objects[ob]->object.c_str();
        dataset_config << '\n';
    }

    // output colors
    dataset_config << "colors\n";
    for (unsigned int ob = 0; ob < objects.size(); ++ob) {
        dataset_config << objects[ob]->color_hex.c_str();
        dataset_config << '\n';
    }

    // close the file
    dataset_config.close();
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	// Parse program arguments
	if (!ParseArgs(argc, argv)) exit(-1);

    // read in all of the gsv scenes
    for (unsigned int ir = 0; ir < run_names.size(); ++ir) {
        char run_name[4096];
        sprintf(run_name, "gsv_data/laser_points/%s/%s.gsv", run_names[ir], run_names[ir]);
        runs.push_back(ReadScene(run_name)->Run(0));
    }
    
    // create the .set file for all of the classes
    if (trainval) CreateSetFile();
    
    // create the training images for each file
    for (unsigned int ir = 0; ir < runs.size(); ++ir) {
        if (DA_images) CreateTraining(runs[ir], "DA");
        if (SA_images) CreateTraining(runs[ir], "SA");
    }

    // save the .set file for all of the classes
    if (trainval) SaveSetFile();

	// Return success 
	return 0;
}