/** \file include/mpi_dispatcher/mpi_skel.hpp
** \brief Declares mpi_skel - a structure to simplify master-slave calculations
*/

#ifndef __INCLUDE_MPISKEL_H
#define __INCLUDE_MPISKEL_H

#include <boost/mpi.hpp>
#include <boost/local_function.hpp>
#include <boost/serialization/vector.hpp>

//#include <type_traits>
#include "mpi_dispatcher.hpp"

#include <memory>

namespace pMPI {

template <typename PartType>
struct ComputeWrap {
    PartType *x;
    int complexity;
    ComputeWrap(PartType &y, int complexity = 1):x(&y),complexity(complexity){};
    void run(){x->compute();};
    ComputeWrap(){};
};

template <typename PartType>
struct PrepareWrap {
    PartType *x;
    int complexity;
    PrepareWrap(PartType &y, int complexity = 1):x(&y),complexity(complexity){};
    void run(){x->prepare();};
    PrepareWrap(){};
};


template <typename WrapType>
struct mpi_skel {
    std::vector<WrapType> parts;
    std::map<pMPI::JobId, pMPI::WorkerId> run(const boost::mpi::communicator& comm, bool VerboseOutput = true);
};

/// Master-slave task schedule, associates jobs with workers and executes mpi_dispatcher to run the simulation
template <typename WrapType>
std::map<pMPI::JobId, pMPI::WorkerId> mpi_skel<WrapType>::run(const boost::mpi::communicator& comm, bool VerboseOutput)
{
    int rank = comm.rank();
    int comm_size = comm.size();
    comm.barrier();
    if (rank==0) { std::cout << "Calculating " << parts.size() << " jobs using " << comm_size << " procs." << std::endl; };

    size_t ROOT = 0;
    std::unique_ptr<pMPI::MPIMaster> disp;

    if (comm.rank() == ROOT) {
        // prepare one Master on a root process for distributing parts.size() jobs
        std::vector<pMPI::JobId> job_order(parts.size());
        for (size_t i=0; i<job_order.size(); i++) job_order[i] = i;
        int BOOST_LOCAL_FUNCTION_TPL(bind this_, std::size_t l, std::size_t r) {
            return (this_->parts[l].complexity > this_->parts[r].complexity); } BOOST_LOCAL_FUNCTION_NAME_TPL(comp1)
        std::sort(job_order.begin(), job_order.end(), comp1);
        disp.reset(new pMPI::MPIMaster(comm,job_order,true));
    };

    comm.barrier();

    // Start calculating data
    for (pMPI::MPIWorker worker(comm,ROOT);!worker.is_finished();) {
        if (rank == ROOT) disp->order();
        worker.receive_order();
        if (worker.is_working()) { // for a specific worker
            JobId p = worker.current_job();
            if (VerboseOutput) std::cout << "["<<p+1<<"/"<<parts.size()<< "] P" << comm.rank()
                                         << " : part " << p << " [" << parts[p].complexity << "] run;" << std::endl;
            parts[p].run();
            worker.report_job_done();
        };
        if (rank == ROOT) disp->check_workers(); // check if there are free workers
    };
    // at this moment all communication is finished
    //comm.barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    // Now spread the information, who did what.
	if (VerboseOutput && rank==ROOT) std::cout << "done." << std::endl;
    comm.barrier();
    std::map<pMPI::JobId, pMPI::WorkerId> job_map;
    if (rank == ROOT) {
        job_map = disp -> DispatchMap;
        std::vector<pMPI::JobId> jobs(job_map.size());
        std::vector<pMPI::WorkerId> workers(job_map.size());

        std::map<pMPI::JobId, pMPI::WorkerId>::const_iterator it = job_map.begin();
        for (int i=0; i<workers.size(); i++) {
            jobs[i] = it->first;
            workers[i] = it->second;
            ++it;
        }
        boost::mpi::broadcast(comm, jobs, ROOT);
        boost::mpi::broadcast(comm, workers, ROOT);
    }
    else {
        std::vector<pMPI::JobId> jobs(parts.size());
        boost::mpi::broadcast(comm, jobs, ROOT);
        std::vector<pMPI::WorkerId> workers(parts.size());
        boost::mpi::broadcast(comm, workers, ROOT);
        for (size_t i=0; i<jobs.size(); i++) job_map[jobs[i]] = workers[i];
    }
    return job_map;
}

}; // end of namespace MPI

#endif
