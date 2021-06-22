/** \file include/mpi_dispatcher/mpi_dispatcher.hpp
** \brief Implementation of a master-slave computation using unblocked MPI Communication that allows to utilize
** master as a computation node
*/
#ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_H
#define POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_H

#include <mpi.h>

#include <stack>
#include <vector>
#include <map>

namespace pMPI {

enum WorkerTag : int { Pending, Work, Finish }; // tags for MPI communication
typedef int JobId;
typedef int WorkerId;

struct MPIWorker
{
    MPI_Comm Comm;
    const WorkerId id;
    const WorkerId boss;

    MPIWorker(const MPI_Comm &comm, WorkerId boss);
    void receive_order();
    void report_job_done();
    bool is_finished();
    bool is_working();

    JobId current_job() const { return current_job_; };

protected:
    JobId current_job_;
    MPI_Request req;
    WorkerTag Status;
};

struct MPIMaster
{
    MPI_Comm Comm;
    size_t Ntasks, Nprocs;

    std::stack<JobId> JobStack;
    std::stack<WorkerId> WorkerStack;

    std::map<JobId, WorkerId> DispatchMap;
    std::vector<JobId> task_numbers;

    std::vector<WorkerId> worker_pool;
    std::map<size_t, WorkerId> WorkerIndices;

    std::vector<MPI_Request> wait_statuses;
    std::vector<bool> workers_finish;

    MPIMaster(const MPI_Comm &comm, std::vector<WorkerId> worker_pool, std::vector<JobId> task_numbers);
    MPIMaster(const MPI_Comm &comm, std::vector<JobId> task_numbers, bool include_boss = true);
    MPIMaster(const MPI_Comm &comm, size_t ntasks, bool include_boss = true);

    void order_worker(WorkerId worker_id, JobId job);
    void order();
    void check_workers();
    bool is_finished() const;

private:
    void fill_stack_();
};

} // namespace pMPI

#endif // #ifndef POMEROL_INCLUDE_MPI_DISPATCHER_MPI_DISPATCHER_H
