/** \file include/mpi_dispatcher/mpi_dispatcher.hpp
** \brief Implementation of a master-slave computation using unbloced MPI Communication that allows to utilize
** master as a computation node
*/


#pragma once

#include <boost/mpi.hpp>
#include <stack>
#include <vector>
#include <map>

namespace pMPI {

enum WorkerTag { Pending, Work, Finish }; // tags for MPI communication
typedef int JobId;
typedef int WorkerId;

struct MPIWorker
{
    boost::mpi::communicator Comm;
    const WorkerId id;
    const WorkerId boss;

    MPIWorker(const boost::mpi::communicator &comm, WorkerId boss);
    void receive_order();
    void report_job_done();
    bool is_finished();
    bool is_working();


    JobId current_job() { return current_job_; };

protected:
    JobId current_job_;
    boost::mpi::request req;
    WorkerTag Status;
};

struct MPIMaster
{
    boost::mpi::communicator Comm;
    size_t Ntasks, Nprocs;

    std::stack<JobId> JobStack;
    std::stack<WorkerId> WorkerStack;

    std::map<JobId, WorkerId> DispatchMap;
    std::vector<JobId> task_numbers;

    std::vector<WorkerId> worker_pool;
    std::map<size_t, WorkerId> WorkerIndices;

    std::vector<boost::mpi::request> wait_statuses;
    std::vector<bool> workers_finish;

    MPIMaster(const boost::mpi::communicator &comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers );
    MPIMaster(const boost::mpi::communicator &comm, std::vector<JobId> task_numbers, bool include_boss = true );
    MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, bool include_boss = true );

    void swap(MPIMaster &x);
    void order_worker(WorkerId worker_id, JobId job);
    void order();
    void check_workers();
    bool is_finished() const;
    void fill_stack_();
};

} // end of namespace MPI


