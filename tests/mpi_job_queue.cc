#include "config.h"
#include <iostream>

#if defined(USING_MPI)

#include "mpi.h"
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

// ----------------------------------------------------------------------------

struct JobResults
{
    int iteration;
    double position[3];
    double dipole[3];
    double duration;
};

double myrandom()
{
    return rand() / (1.0 + RAND_MAX);
}

int randint(int a, int b)
{
    return a + (int)ceil((b - a) * myrandom());
}

// ----------------------------------------------------------------------------

#define READY_FOR_JOB   1
#define NEW_JOB         2
#define JOB_UPDATE      3
#define NO_MORE_JOBS    4

void send_int(int msg, int dest, int tag)
{
    MPI::COMM_WORLD.Send(&msg, 1, MPI::INT, dest, tag);
}

void receive_int(int &msg, int& source, int& tag)
{
    MPI::Status status;
    MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
    source = status.Get_source();
    tag = status.Get_tag();
}

void send_job_results(JobResults& data, int dest)
{
    int ibuf[1];
    ibuf[0] = data.iteration;
    MPI::COMM_WORLD.Send(ibuf, 1, MPI::INT, dest, JOB_UPDATE);

    double dbuf[7];
    dbuf[0] = data.position[0];
    dbuf[1] = data.position[1];
    dbuf[2] = data.position[2];
    dbuf[3] = data.dipole[0];
    dbuf[4] = data.dipole[1];
    dbuf[5] = data.dipole[2];
    dbuf[6] = data.duration;
    MPI::COMM_WORLD.Send(dbuf, 7, MPI::DOUBLE, dest, JOB_UPDATE);
}

void receive_job_results(JobResults& data, int source)
{
    MPI::Status status;

    int ibuf[1];
    MPI::COMM_WORLD.Recv(ibuf, 1, MPI::INT, source, JOB_UPDATE);
    data.iteration = ibuf[0];

    double dbuf[7];
    MPI::COMM_WORLD.Recv(dbuf, 7, MPI::DOUBLE, source, JOB_UPDATE);
    data.position[0] = dbuf[0];
    data.position[1] = dbuf[1];
    data.position[2] = dbuf[2];
    data.dipole[0] = dbuf[3];
    data.dipole[1] = dbuf[4];
    data.dipole[2] = dbuf[5];
    data.duration = dbuf[6];
}

// ----------------------------------------------------------------------------

void run_job(int id)
{
    JobResults data;
    int master = 0;

    // XXX: for now, simulate a job callback
    int n_iterations = randint(10, 50);
    for (int i = 0; i < n_iterations; i++)
    {
        sleep(randint(1,3));

        send_int(id, master, JOB_UPDATE);

        data.iteration = i;

        data.position[0] = myrandom();
        data.position[1] = myrandom();
        data.position[2] = myrandom();

        data.dipole[0] = myrandom();
        data.dipole[1] = myrandom();
        data.dipole[2] = myrandom();

        data.duration = myrandom();

        send_job_results(data, master);
    }
}

void run_worker()
{
    using namespace std;

    // rank of master process
    int master = 0;

    // the worker process rank
    int rank = MPI::COMM_WORLD.Get_rank();

    cout << rank << ": Starting worker" << endl;

    int msg, tag;
    MPI::Status status;

    // seed random number generator
    srand(time(NULL) + rank);

    while (true)
    {
        // ready for another one
        send_int(0, master, READY_FOR_JOB);

        // check for a job
        MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, master, MPI::ANY_TAG, status);
        tag = status.Get_tag();
        if (tag == NEW_JOB)
        {
            int id = msg;
            run_job(id);
        }
        else if (tag == NO_MORE_JOBS)
        {
            break;
        }
    }

    return;
}

void run_controller()
{
    using namespace std;

    int rank = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();
    const int num_workers = nprocs-1;

    assert(rank == 0);
    cout << rank << ": Starting controller" << endl;

    JobResults data;
    int msg, src, tag;

    // build the jobs
    const int num_jobs = 10;
    typedef std::vector<int> JobCollection;
    JobCollection jobs;
    for (int i = 0; i < num_jobs; i++)
        jobs.push_back(i);

    // count the idle workers (once all jobs are done)
    int finished = 0;

    // dispatch the jobs
    while (true)
    {
        if (finished == num_workers)
        {
            cout << rank << ": All jobs finished!" << endl;
            break;
        }

        receive_int(msg, src, tag);

        if (tag == READY_FOR_JOB)
        {
            if (jobs.empty())
            {
                cout << src << ": Done" << endl;
                send_int(0, src, NO_MORE_JOBS);
                finished++;
                continue;
            }

            int job_id = jobs.back();
            jobs.pop_back();

            cout << rank << ": Sending job " << job_id << " to process " << src << endl;
            send_int(job_id, src, NEW_JOB);
        }
        else if (tag == JOB_UPDATE)
        {
            int id = msg;
            receive_job_results(data, src);

            cout << src << ": job " << id
                 << ", iteration " << data.iteration
                 << ", x = [" << data.position[0] << ", " << data.position[1] << ", " << data.position[2] << "]"
                 << ", p = [" << data.dipole[0] << ", " << data.dipole[1] << ", " << data.dipole[2] << "]"
                 << ", T = " << data.duration << " secs"
                 << endl;
        }
    }
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == 0)
        run_controller();
    else
        run_worker();
    MPI::Finalize();
    return 0;
}

// ----------------------------------------------------------------------------
#else

int main(void)
{
    std::cerr << "MPI not available" << std::endl;
    return 1;
}

#endif // USING_MPI
