#include "MPIDispatcher.h"

namespace Pomerol { 

namespace MPI {

//
// Worker
//

bool MPIWorker::is_finished()
{
if (Status == WorkerTag::Finish) return true; // guard to avoid double testing
Status = (FinishReq.test()?WorkerTag::Finish:Status);
WorkReq = boost::mpi::request();
if (Status == WorkerTag::Finish) return true;
return false;
}

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


MPIMaster::MPIMaster(const boost::mpi::communicator &comm, size_t ntasks, bool include_boss):
    Comm(comm),Ntasks(ntasks),
    Nprocs(comm.size()-!include_boss),id(Comm.rank()),
    task_numbers(Ntasks),
    wait_statuses(Nprocs)
{
    if (!Nprocs) throw (std::logic_error("No workers to evaluate"));
    for (size_t i=0; i<Ntasks; i++) { 
        task_numbers[i] = i; 
        JobStack.push(task_numbers[i]); 
    };
    for (size_t p=0; p<comm.size(); p++) { 
        if (include_boss || id != p) {
            worker_pool.push_back(p); 
            WorkerIndices[p] = p; 
            WorkerStack.push(p); 
            };
        };
};



void MPIMaster::order_worker(ProcId worker, JobId job)
{
    DEBUG("Occupying worker " << worker << " with job = " << job);
    Comm.isend(worker,int(WorkerTag::Work),job);
    wait_statuses[WorkerIndices[worker]] = Comm.irecv(worker,int(WorkerTag::Pending));
};

void MPIMaster::order()
{
    while (!WorkerStack.empty()) { 
        auto worker = WorkerStack.top(); 
        auto job = JobStack.top(); 
        order_worker(worker,job); 
        WorkerStack.pop(); 
        JobStack.pop(); 
    }; 
};

void MPIMaster::update()
{
    INFO("!");
    if (!JobStack.empty()) {
        for (size_t i=0; i<Nprocs; i++) {
            if (wait_statuses[i].test()) { 
                WorkerStack.push(worker_pool[i]);
                wait_statuses[i] = boost::mpi::request();
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

