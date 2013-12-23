#include "MPIDispatcher.h"

namespace Pomerol { 

namespace MPI {

//
// Worker
//
MPIWorker::MPIWorker(const boost::mpi::communicator &comm, ProcId boss):
        Comm(comm),
        id(Comm.rank()),
        boss(boss),
        WorkReq(Comm.irecv(boss, int(WorkerTag::Work), current_job)),
        FinishReq(Comm.irecv(boss, int(WorkerTag::Finish))),
        Status(WorkerTag::Pending)
    {};



bool MPIWorker::is_finished()
{
return (Status == WorkerTag::Finish);
}

bool MPIWorker::is_working()
{
return (Status == WorkerTag::Work);
}

void MPIWorker::receive_order() 
{
    if (Status == WorkerTag::Pending && WorkReq.test()) { Status = WorkerTag::Work; }; 
    if (Status == WorkerTag::Pending && FinishReq.test()) { Status = WorkerTag::Finish; };
}

void MPIWorker::report_job_done()
{ 
    Comm.isend(boss,int(WorkerTag::Pending)); 
    Status = WorkerTag::Pending; 
    WorkReq = Comm.irecv(boss, int(WorkerTag::Work), current_job);
};
//
// Master
//

MPIMaster::MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, std::vector<ProcId> worker_pool, std::vector<JobId> task_numbers ):
    Comm(comm), Ntasks(ntasks),
    Nprocs(worker_pool.size()), id(Comm.rank()),
    task_numbers(task_numbers), worker_pool(worker_pool), 
    wait_statuses(Nprocs)
{
    for (size_t i=0; i<Ntasks; i++) { JobStack.push(task_numbers[i]); };
    for (size_t p=0; p<Nprocs; p++) { 
        WorkerIndices[worker_pool[p]] = p; 
        WorkerStack.push(worker_pool[p]); 
        };
};

inline std::vector<ProcId> _autorange_workers(const boost::mpi::communicator &comm, bool include_boss)
{
    std::vector<ProcId> out;
    size_t Nprocs(comm.size()-!include_boss);
    if (!Nprocs) throw (std::logic_error("No workers to evaluate"));
    for (size_t p=0; p<comm.size(); p++) { 
        if (include_boss || comm.rank() != p) {
            out.push_back(p); 
            };
        };
    return out;
}

inline std::vector<JobId> _autorange_tasks(size_t ntasks)
{
    std::vector<JobId> out(ntasks);
    for (size_t i=0; i<ntasks; i++) { 
        out[i] = i; 
    };
    return out;
}

MPIMaster::MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, bool include_boss):
    MPIMaster(comm,ntasks,_autorange_workers(comm,include_boss), _autorange_tasks(ntasks))
{};


void MPIMaster::order_worker(ProcId worker, JobId job)
{
    DEBUG("Occupying worker " << worker << " with job = " << job);
    Comm.isend(worker,int(WorkerTag::Work),job);
    wait_statuses[WorkerIndices[worker]] = Comm.irecv(worker,int(WorkerTag::Pending));
};

void MPIMaster::order()
{
    while (!WorkerStack.empty() && !JobStack.empty()) { 
        auto worker = WorkerStack.top(); 
        auto job = JobStack.top(); 
        order_worker(worker,job); 
        WorkerStack.pop(); 
        JobStack.pop(); 
    }; 
};

void MPIMaster::update()
{
    if (!JobStack.empty()) {
        for (size_t i=0; i<Nprocs && !JobStack.empty(); i++) {
            if (wait_statuses[i].test()) { 
                WorkerStack.push(worker_pool[i]);
        //        wait_statuses[i] = boost::mpi::request();
                }; 
        };
    }
    else {
        for (size_t i=0; i<Nprocs; i++) {
            Comm.isend(worker_pool[i],int(WorkerTag::Finish));
        };
    };
}

} // end of namespace MPI
} // end of namespace Pomerol

