// Takes in a script file via command line and outputs new scripts for each
// processor. Those scripts are also run.

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

#if defined(USING_MPI)
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes
    
    ifstream scriptfile (argv[1]); // format: mpi_jobarray blah.sh??

    int N = 0;
    string line;
    // Count the number of lines.
    while (getline(scriptfile, line))
    {
	N++;
    }
    scriptfile.clear();
    scriptfile.seekg(0, ios::beg);
    
    
    // Probably a better way to do this, but eh.
    int nlines = N-1; // Don't count the first line.
    vector<int> nvector (size, nlines / size); // initial assignment
    int remainder = nlines - ((nlines / size) * size);
    // Distribute the remaining lines.
    for (int i = 0; i < remainder; i++)
    {
	nvector[i]++;
    }
    // Prepare indices for each processor.
    int begin = 1;
    for (int i = 0; i < rank; i++)
    {
	begin += nvector[i];
    }
    int line_number = 0;
    vector<string> exec_lines;
    while (getline(scriptfile, line))
    { 
	if ((line_number >= begin) && (line_number < begin + nvector[rank]))
	{
	    exec_lines.push_back(line);
	}
	line_number++;
    }

    stringstream ss;
    ss << "script" << rank << ".sh";
    string filename = ss.str();
    ofstream outfile (filename.c_str());

    // Write the new script.
    outfile << "#!/bin/bash -x" << endl;
    
    for (int i = 0; i < nvector[rank]; i++)
    {
	outfile << exec_lines[i] << endl;
    }

    outfile.close();
    scriptfile.close();

    // Run the script.
    string command = "chmod +x " + filename;
    system(command.c_str());
    system(("./"+filename).c_str());

    MPI_Finalize();
    return 0;
}

#else

int main(void)
{
    fprintf(stderr, "MPI not available\n");
    return 1;
}

#endif // USING_MPI
