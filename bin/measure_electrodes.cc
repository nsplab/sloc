/*
 * Given a potentials data file, and an electrodes file, write out a electrode measurements file
 * that we can use in bem_inverse_solve. This is a fake measurement (we're only adding random noise)
 * for testing purposes.
 */

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <sloc/random.h>

using namespace std;

// ----------------------------------------------------------------------------

double rms(vector<double>& x)
{
    const int n = x.size();
    double sum2 = 0;
    for (int i = 0; i < n; i++)
        sum2 += x[i] * x[i];
    return sqrt(sum2 / n);
}

// ----------------------------------------------------------------------------

void usage(const char *pgm)
{
    cerr << "Usage: " << pgm
         << " -p datfile -e electrodes [-n SNR] -o outfile"
         << endl;
    exit(1);
}

void process_args(int argc, char *argv[], string& datfile, string& ecfile, string& outfile, double& SNR)
{
    char *pgm = argv[0];

    if (argc == 1)
        usage(pgm);

    ecfile = "";
    datfile = "";
    outfile = "";
    SNR = numeric_limits<double>::infinity();

    while (argc > 1)
    {
        if (strncmp(argv[1], "-p", 3) == 0)
        {
            // read the next argument (expecting datfile)
            argc--; argv++;
            if (argc > 1)
                datfile = argv[1];
            else
                break;
        }
        else if (strncmp(argv[1], "-e", 3) == 0)
        {
            // read the next argument (expecting electrodes file)
            argc--; argv++;
            if (argc > 1)
                ecfile = argv[1];
            else
                break;
        }
        else if (strncmp(argv[1], "-n", 3) == 0)
        {
            // read the next argument (expecting snr noise level)
            argc--; argv++;
            if (argc > 1)
                try {
                    SNR = boost::lexical_cast<double>(argv[1]);
                } catch(boost::bad_lexical_cast &) {
                    cerr << "Invalid SNR '" << argv[1] << "'" << endl;
                    usage(pgm);
                }
            else
                break;
        }
        else if (strncmp(argv[1], "-o", 3) == 0)
        {
            // read the next argument (expecting outfile)
            argc--; argv++;
            if (argc > 1)
                outfile = argv[1];
            else
                break;
        }

        // advance to the next argument
        argc--; argv++;
    }

    if (datfile.empty())
    {
        cerr << "Missing data file!" << endl;
        usage(pgm);
    }

    if (ecfile.empty())
    {
        cerr << "Missing electrodes file!" << endl;
        usage(pgm);
    }

    if (outfile.empty())
    {
        cerr << "Missing output file!" << endl;
        usage(pgm);
    }
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    ifstream in;
    ofstream out;
    string datfile, ecfile, outfile;
    double SNR;
    int i;

    // process cli args
    process_args(argc, argv, datfile, ecfile, outfile, SNR);

    // read the data file
    vector<double> phi;
    cout << "Reading potentials '" << datfile << "'" << endl;
    in.open(datfile.c_str());
    if (!in.is_open())
    {
        cerr << "Could not open potentials data file '" << datfile << "'" << endl;
        usage(argv[0]);
    }
    while (!in.eof())
    {
        string line;

        getline(in, line);
        if (line.empty())
            break;

        try {
            phi.push_back(boost::lexical_cast<double>(line));
        } catch (boost::bad_lexical_cast &) {
            cerr << "Bad line in potentials data file: " << line << endl;
            exit(2);
        }
    }
    in.close();

    // read the electrodes
    vector<unsigned int> ec;
    cout << "Reading electrodes '" << ecfile << "'" << endl;
    in.open(ecfile.c_str());
    if (!in.is_open())
    {
        cerr << "Could not open electrodes file '" << ecfile << "'" << endl;
        usage(argv[0]);
    }
    while (!in.eof())
    {
        string line;

        getline(in, line);
        if (line.empty())
            break;

        try {
            ec.push_back(boost::lexical_cast<unsigned int>(line));
        } catch (boost::bad_lexical_cast &) {
            cerr << "Bad line in electrodes file: " << line << endl;
            exit(2);
        }
    }
    in.close();

    // XXX: shouldn't the noise level depend on the layer? think about this...
    // for now, just assume that the noise is applied equally to all nodes.

    // generate the noise (save it to a file?)
    const double phi_rms = rms(phi);
    const double snr_factor = 1/sqrt(SNR);
    const double sigma = phi_rms * snr_factor;
    cout << "Calculated phi_rms = " << phi_rms << endl;
    cout << "Using snr = " << SNR << ", snr_factor = " << snr_factor << endl;
    cout << "Generating noise with sigma " << sigma << endl;
    const int num_electrodes = ec.size();
    double noise[num_electrodes];
    sloc::seed_random_generator();
    for (i = 0; i < num_electrodes; i++)
        noise[i] = sloc::normal(0, sigma);

    // apply the noise
    double em[num_electrodes];
    for (i = 0; i < num_electrodes; i++)
        em[i] = phi[ec[i]] + noise[i];

    // write out the electrode measurements file
    out.open(outfile.c_str());
    out << num_electrodes << endl;
    cout << "Processing " << num_electrodes << " electrode measurements to '" << outfile << "'" << endl;
    for (i = 0; i < num_electrodes; i++)
        out << ec[i] << " " << em[i] << endl;
    out.close();

    // done!
    return 0;
}
