// Source file for the simple dijkstra algorithm

// Include files 

#include "Neuron/Neuron.h"

// program arguments
static int print_verbose = 0;
static int print_debug = 0;
static int supervoxels = 0;
static char *affinity_filename = NULL;
static char *image_filename = NULL;
static char *machine_filename = NULL;

// global variables
static NeuronData *nd;

////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void) {
	// Start statistics
	RNTime start_time;
	start_time.Read();

	// allocate new neuron data
	nd = new NeuronData();

	// read in the various files
	if (!nd->ReadAffinityMetaFile(affinity_filename)) { fprintf(stderr, "Failed to read %s.\n", affinity_filename); return 0; }
	if (image_filename && !nd->ReadImageMetaFile(image_filename)) { fprintf(stderr, "Failed to read %s.\n", image_filename); return 0; }
	if (machine_filename && !nd->ReadMachineMetaFile(machine_filename)) { fprintf(stderr, "Failed to read %s.\n", machine_filename); return 0; }

	// print statistics
	if (print_verbose) {
		printf("Read grid...\n");
		printf("  Time = %.2f seconds\n", start_time.Elapsed());
		fflush(stdout);
	}

	// return success
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv) {
    // Parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) print_debug = 1;
			else if (!strcmp(*argv, "-affinity")) { argv++; argc--; affinity_filename = *argv; }
			else if (!strcmp(*argv, "-image")) { argv++; argc--; image_filename = *argv; }
			else if (!strcmp(*argv, "-machine")) { argv++; argc--; machine_filename = *argv; }
			else if (!strcmp(*argv, "-supervoxels")) { supervoxels = 1; }
			else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
        }
        argv++; argc--;
    }

	if (!affinity_filename) { fprintf(stderr, "Need to supply the affinities filename.\n"); return 0; }

	if (print_verbose) {
		if (!image_filename) {
			printf("No image filename supplied (OK to continue).\n");
		}
		if (!machine_filename) {
			printf("Cannot use supervoxels.\n");
			supervoxels = 0;
		}
	}

    // Return OK status 
    return 1;
}

////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    // Parse program arguments
    if (!ParseArgs(argc, argv)) exit(-1);

	// read in the neuron data
	if (!ReadData()) exit(-1);

	RNArray<NeuronVoxel *> nv = RNArray<NeuronVoxel *>();

    // Return success
    return 0;
}
