from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
data = comm.gather(f"hello from rank {rank}/{size}", root=0)
if rank == 0:
    print(data)
