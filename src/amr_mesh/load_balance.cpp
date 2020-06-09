/*!
  @file load_balance.cpp
  @details Implementatino of load_balance.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "load_balance.H"
#include "EBEllipticLoadBalance.H"

void load_balance::balance_volume(Vector<int>& a_procs, const Vector<Box>& a_boxes){
  LoadBalance(a_procs, a_boxes);
}

void load_balance::balance_elliptic(Vector<int>&                       a_procs,
				    const Vector<Box>&                 a_boxes,
				    const RefCountedPtr<EBIndexSpace>& a_ebis,
				    const ProblemDomain&               a_domain,
				    const bool                         a_verbose){

  EBEllipticLoadBalance(a_procs, a_boxes, a_domain, a_verbose, a_ebis);
}

void load_balance::balance_multifluid(Vector<int>&               a_procs,
				      const Vector<Box>&         a_boxes,
				      const RefCountedPtr<mfis>& a_mfis,
				      const ProblemDomain&       a_domain,
				      const bool                 a_verbose){

  MayDay::Abort("load_balance::balance_multifluid - not implemented (yet)");
  
  // Load balance the conventional way
  Vector<Box> boxes = a_boxes;
  LoadBalance(a_procs, a_boxes);

  
  MayDay::Abort("load_balance::balance_multifluid - not implemented");
}

void load_balance::get_mfpoisson_loads(Vector<unsigned long long>& a_loads,
				       Vector<Box>&                a_boxes,
				       RefCountedPtr<mfis>&        a_mfis,
				       const DisjointBoxLayout&    a_dbl,
				       const ProblemDomain&        a_domain){
  
}

void load_balance::get_cdr_loads(Vector<unsigned long long>& a_loads,
				 Vector<Box>&                a_boxes,
				 RefCountedPtr<mfis>&        a_mfis,
				 const DisjointBoxLayout&    a_dbl,
				 const ProblemDomain&        a_domain){
  
}

void load_balance::get_rte_loads(Vector<unsigned long long>& a_loads,
				 Vector<Box>&                a_boxes,
				 RefCountedPtr<mfis>&        a_mfis,
				 const DisjointBoxLayout&    a_dbl,
				 const ProblemDomain&        a_domain){
  
}

void load_balance::gather_boxes(Vector<Box>& a_boxes){

  // TLDR: This code does a gather operation on the loads and boxes. They are gather globally in this way:
  //       (BoxesForMPIRank=0, BoxesForMPIRank=1, BoxesForMPIRank=2, ....)

  // 1. Linearize local boxes
  int send_size    = 2*CH_SPACEDIM;                 // Message size for one box
  int send_count   = a_boxes.size()*send_size;      // Number of elements sent from this rank
  int* send_buffer = new int[send_count];           // Send buffer for this MPI rank
  int* send_buf2   = send_buffer;                   // Backup address. Going to monkey with pointer increments on send buffer

  // Linearize a_boxes onto send_buffer
  for (int i = 0; i < a_boxes.size(); i++, send_buffer+=send_size){
    const Box& b = a_boxes[i];
    D_TERM6(send_buffer[0] =b.smallEnd(0); send_buffer[1] =b.bigEnd(0);,
	    send_buffer[2] =b.smallEnd(1); send_buffer[3] =b.bigEnd(1);,
	    send_buffer[4] =b.smallEnd(2); send_buffer[5] =b.bigEnd(2);,
	    send_buffer[6] =b.smallEnd(3); send_buffer[7] =b.bigEnd(3);,
	    send_buffer[8] =b.smallEnd(4); send_buffer[9] =b.bigEnd(4);,
	    send_buffer[10]=b.smallEnd(5); send_buffer[11]=b.bigEnd(5););
  }
  send_buffer = send_buf2; // Revert point to start of array


  // 2. Get the number of elements sent from each rank
  int* send_counts = new int[numProc()];
  MPI_Allgather(&send_count, 1, MPI_INT, send_counts, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0] = 0;
  for (int i = 0; i < numProc()-1; i++){
    offsets[i+1] = offsets[i] + send_counts[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++){
    total_count += send_counts[i];
  }
  int* recv_buffer = new int[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, send_count, MPI_INT, recv_buffer, send_counts, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_boxes.resize(0);
  int* recv_buf2 = recv_buffer; // Going to monkey with pointer increments again
  for (int i = 0; i < total_count/send_size; i++, recv_buffer+=send_size){
    IntVect lo, hi;
    D_TERM6(lo[0] = recv_buffer[0];  hi[0] = recv_buffer[1];,
	    lo[1] = recv_buffer[2];  hi[1] = recv_buffer[3];,
	    lo[2] = recv_buffer[4];  hi[2] = recv_buffer[5];,
	    lo[3] = recv_buffer[6];  hi[3] = recv_buffer[7];,
	    lo[4] = recv_buffer[8];  hi[4] = recv_buffer[9];,
	    lo[5] = recv_buffer[10]; hi[5] = recv_buffer[11];);

    a_boxes.push_back(Box(lo, hi));
  }
  recv_buffer = recv_buf2;

  delete recv_buffer;
  delete send_buffer;

}

void load_balance::gather_loads(Vector<Real>& a_loads){

  // TLDR: This code does a gather operation on the loads. They are gather globally in this way:
  //       (LoadsForRank=0, LoadsForRank=1, LoadsForRank2=2, ....)

  // 1. Linearize local boxes onto appropriate memory structure
  int send_size     = 1;                             // Message size for one element. Just one Real. 
  int mySendCount    = a_loads.size()*send_size;      // Number of elements sent from this rank
  Real* send_buffer = new Real[mySendCount];         // Send buffer for this MPI rank
  for (int i = 0; i < a_loads.size(); i++){
    send_buffer[i] = a_loads[i];
  }

  // 2. Get the number of elements sent from each rank
  int* allSendCount = new int[numProc()];
  MPI_Allgather(&mySendCount, 1, MPI_INT, allSendCount, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0] = 0;
  for (int i = 0; i < numProc()-1; i++){
    offsets[i+1] = offsets[i] + allSendCount[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++){
    total_count += allSendCount[i];
  }
  Real* recv_buffer = new Real[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, mySendCount, MPI_CH_REAL, recv_buffer, allSendCount, offsets, MPI_CH_REAL, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_loads.resize(total_count);
  for (int i = 0; i < total_count; i++){
    a_loads[i] = recv_buffer[i];
  }

  delete recv_buffer;
  delete send_buffer;
}

void load_balance::gather_loads(Vector<int>& a_loads){

  // TLDR: This code does a gather operation on the loads. They are gather globally in this way:
  //       (LoadsForRank=0, LoadsForRank=1, LoadsForRank2=2, ....)

  // 1. Linearize local boxes onto appropriate memory structure
  int send_size     = 1;                             // Message size for one element. Just one int. 
  int mySendCount    = a_loads.size()*send_size;      // Number of elements sent from this rank
  int* send_buffer = new int[mySendCount];         // Send buffer for this MPI rank
  for (int i = 0; i < a_loads.size(); i++){
    send_buffer[i] = a_loads[i];
  }

  // 2. Get the number of elements sent from each rank
  int* allSendCount = new int[numProc()];
  MPI_Allgather(&mySendCount, 1, MPI_INT, allSendCount, 1, MPI_INT, Chombo_MPI::comm);

  // 3. Compute offsets
  int* offsets = new int[numProc()];
  offsets[0] = 0;
  for (int i = 0; i < numProc()-1; i++){
    offsets[i+1] = offsets[i] + allSendCount[i];
  }

  // 4. Allocate storage for total buffer size
  int total_count = 0;
  for (int i = 0; i < numProc(); i++){
    total_count += allSendCount[i];
  }
  int* recv_buffer = new int[total_count];

  // 5. MPI send
  MPI_Allgatherv(send_buffer, mySendCount, MPI_INT, recv_buffer, allSendCount, offsets, MPI_INT, Chombo_MPI::comm);

  // 6. Delinearize buffer, make it into boxes
  a_loads.resize(total_count);
  for (int i = 0; i < total_count; i++){
    a_loads[i] = recv_buffer[i];
  }

  delete recv_buffer;
  delete send_buffer;
}

void load_balance::gather_boxes_and_loads(Vector<Box>& a_boxes, Vector<int>& a_loads){
  load_balance::gather_boxes(a_boxes);
  load_balance::gather_loads(a_loads);
}

void load_balance::load_balance_boxes(Vector<int>&       a_procs,
				      const Vector<int>& a_loads,
				      const Vector<Box>& a_boxes){

  if(a_loads.size()!= a_boxes.size()){
    MayDay::Abort("load_balance::load_balance_boxes - logic bust, loads.size() != boxes.size()");
  }
  
  Vector<long int> loads(a_loads.size());
  for (int i = 0; i < loads.size(); i++){
    loads[i] = a_loads[i];
  }

  LoadBalance(a_procs, loads, a_boxes);
}
