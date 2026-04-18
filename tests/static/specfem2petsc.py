import numpy as np
from mpi4py import MPI 

def get_hash_score(iglob: int, myrank: int) -> int:
    mask_64 = 0xFFFFFFFFFFFFFFFF
    iglob = int(iglob)
    myrank = int(myrank)

    # Use Python integers plus an explicit 64-bit mask so the arithmetic
    # matches uint64 wraparound without NumPy overflow warnings.
    score = ((iglob * 0x9E3779B97F4A7C15) ^ (myrank * 0xBF58476D1CE4E5B9)) & mask_64

    score ^= score >> 30
    score = (score * 0xBF58476D1CE4E5B9) & mask_64

    score ^= score >> 27
    score = (score * 0x94D049BB133111EB) & mask_64

    score ^= score >> 31
    score &= mask_64
    return int(score)

def mpi_print(*args, **kwargs) -> None:
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    size = comm.Get_size()

    for i in range(size):
        if myrank == i:
            print(f"Rank {myrank}: ", *args, **kwargs)
        comm.Barrier() # synchronize before the next rank prints

class Mesh:
    def __init__(self, nspec=2, ngll=4):

        comm = MPI.COMM_WORLD
        self.myrank = comm.Get_rank()
        
        self.nspec = nspec
        self.ngll = ngll

        # connectivity array ibool, shape (nspec, ngll), stores the global node indices for each element, this is just a dummy array for testing purposes, in practice, this should be generated based on the mesh partitioning and connectivity
        self.ibool = np.zeros((nspec, ngll), dtype=int)
        self.ibool[0,:] = [0,1,2,3]
        self.ibool[1,:] = [2,3,4,5]
        self.nglob = np.max(self.ibool) + 1 # total number of global nodes, which is 6 in this case

        # coordinates
        self.x = np.zeros((self.nglob), dtype=float)
        self.y = np.zeros((self.nglob), dtype=float)
        self.z = np.zeros((self.nglob), dtype=float)

        # Two adjacent 2D elements: rank 0 spans x in [-1, 0], rank 1 spans x in [0, 1].
        # Both partitions use the same y coordinates, so nodes 2 and 3 are shared.
        x_coords = np.array([-1.0, 0, -1, 0.0, -1.0, 0.0], dtype=float)
        y_coords = np.array([-5.0,  -5.0, 0.0, 0, 5.0, 5.0], dtype=float)

        for ispec in range(self.nspec):
            for ipt in range(self.ngll):
                iglob = self.ibool[ispec, ipt]
                self.x[iglob] = x_coords[iglob] + self.myrank
                self.y[iglob] = y_coords[iglob]
                self.z[iglob] = 0.0

        # data exchange information for parallel assembly
        self.num_neighbors = 1
        self.neighbor_ranks = np.zeros(self.num_neighbors, dtype=int)
        self.neighbor_ranks[0] = 1 - self.myrank # rank 0's neighbor is rank 1, and vice versa
        self.max_send_size = 3
        self.neighbor_points = np.zeros(self.num_neighbors, dtype=int)
        self.neighbor_points[0] = self.max_send_size # for simplicity, we assume each neighbor requires the same number of points to be sent, in practice, this should be determined by the mesh partitioning and connectivity
        self.send_points = np.zeros((self.num_neighbors, self.max_send_size), dtype=int)
        if self.myrank == 0:
            self.send_points[0,:] = [1,3,5] # rank 0 needs to send data for nodes 1,3,5 to rank 1
        else:
            self.send_points[0,:] = [0,2,4] # rank 1 needs to send data for nodes 0,2,4 to rank 0

        # create a sparse matrix for data exchange
        self.xadj_neigbors = np.zeros((self.num_neighbors + 1), dtype=int)
        for ib in range(self.num_neighbors):
            self.xadj_neigbors[ib + 1] = self.xadj_neigbors[ib] + self.neighbor_points[ib]
        self.neighbor_adj = np.zeros((self.xadj_neigbors[-1]), dtype=int)
        for ib in range(self.num_neighbors):
            istart = self.xadj_neigbors[ib]
            npts = self.neighbor_points[ib]
            for ipt in range(npts):
                self.neighbor_adj[istart + ipt] = self.send_points[ib, ipt]

    def build(self):
        self.owner_ranks = np.zeros(self.nglob, dtype=int)
        self.owner_ranks[:] = self.myrank # for simplicity, we assign all nodes

        # lower rank owns the shared nodes, in practice, this should be determined by the mesh partitioning
        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            npts = self.neighbor_points[ib]
            for ipt in range(npts):
                iglob = self.send_points[ib, ipt]

                # get coordinates
                x0 = self.x[iglob]
                y0 = self.y[iglob]
                z0 = self.z[iglob]

                # get a big integer based on coordinates
                hash_val = int((x0 * 73856093) + (y0 * 19349663) + (z0 * 83492791)) # a common hash function for 3D coordinates, the large prime numbers help to distribute the hash values more uniformly

                score_current = get_hash_score(hash_val, self.owner_ranks[iglob])
                score_neighbor = get_hash_score(hash_val, nrank)
                if score_neighbor > score_current: # if the neighbor's score is better, then the neighbor rank becomes the new owner of this node
                    self.owner_ranks[iglob] = nrank

                # if nrank < self.owner_ranks[iglob]: # if the neighbor rank is less than the current owner rank, then the neighbor rank becomes the new owner of this node
                #     self.owner_ranks[iglob] = nrank 
        self.nglob_owned = np.sum(self.owner_ranks == self.myrank) # count how many nodes are owned by the current rank
        #print(f"Rank {self.myrank}: nglob = {self.nglob}, nglob_owned = {self.nglob_owned}, owner_ranks = {self.owner_ranks}")
        mpi_print(f"nglob = {self.nglob}, nglob_owned = {self.nglob_owned}, owner_ranks = {self.owner_ranks}")

    def _create_l2g_map(self):
        comm = MPI.COMM_WORLD
        global_offset = comm.exscan(self.nglob_owned)
        if self.myrank == 0:
            global_offset = 0

        self.l2g_map = np.zeros(self.nglob, dtype=int)
        owner_idx = 0
        for iglob in range(self.nglob):
            if self.owner_ranks[iglob] == self.myrank:
                self.l2g_map[iglob] = global_offset + owner_idx
                owner_idx += 1
            else:
                self.l2g_map[iglob] = -1 # not owned by this rank

        # data exchange to fill in the l2g_map entries for non-owned nodes, in practice, this should be done using non-blocking MPI communication to overlap with computation
        buf_sd = np.zeros((1, self.max_send_size), dtype=int)
        buf_rv = np.zeros((1, self.max_send_size), dtype=int)
        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            npts = self.neighbor_points[ib]
            for ipt in range(npts):
                iglob = self.send_points[ib, ipt]
                buf_sd[ib, ipt] = self.l2g_map[iglob] # send the l2g_map entry for this node to the neighbor
                comm.Sendrecv(sendbuf=buf_sd, dest=nrank, recvbuf=buf_rv, source=nrank)

        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            npts = self.neighbor_points[ib]
            for ipt in range(npts):
                iglob = self.send_points[ib, ipt]
                if self.owner_ranks[iglob] == nrank : # if this node is not owned by the current rank, then update the l2g_map entry with the received value from the neighbor
                    self.l2g_map[iglob] = buf_rv[ib, ipt]

        # get start/endid 
        self.rstart = global_offset
        self.rend = global_offset + self.nglob_owned

        mpi_print(f"l2g_map = {self.l2g_map}, rstart = {self.rstart}, rend = {self.rend}")

    def _create_node_adjacency(self):
        # 1. First create node2elem adjacency counts safely
        num_elems_per_node = np.zeros((self.nglob), dtype=int)
        for ispec in range(self.nspec):
            for ipt in range(self.ngll):
                iglob = self.ibool[ispec, ipt]
                num_elems_per_node[iglob] += 1

        mpi_print(f"num_elems_per_node = {num_elems_per_node}")
        
        # Build node2elem pointers
        node2elem_ptr = np.zeros((self.nglob + 1), dtype=int)
        for iglob in range(self.nglob):
            node2elem_ptr[iglob + 1] = node2elem_ptr[iglob] + num_elems_per_node[iglob]

        # Populate node2elem data
        node2elem_data = np.zeros((node2elem_ptr[-1]), dtype=int)
        node_elem_count = np.zeros((self.nglob), dtype=int)
        for ispec in range(self.nspec):
            for ipt in range(self.ngll):
                iglob = self.ibool[ispec, ipt]
                idx = int(node2elem_ptr[iglob] + node_elem_count[iglob])
                node2elem_data[idx] = ispec
                node_elem_count[iglob] += 1

        # 2. Create node adjacency pointers (Cumulative sum built on the fly!)
        self.node_adj_ptr = np.zeros((self.nglob + 1), dtype=int)
        node_mask = np.zeros((self.nglob), dtype=int) 
        node_mask[:] = -1 
        
        for iglob in range(self.nglob):
            istart = node2elem_ptr[iglob]
            iend = node2elem_ptr[iglob + 1]
            
            # *** TYPO FIXED HERE *** # Carry over the previous cumulative total before adding new neighbors
            self.node_adj_ptr[iglob + 1] = self.node_adj_ptr[iglob] 
            
            for idx in range(istart, iend):
                ispec = node2elem_data[idx]
                for ipt in range(self.ngll):
                    iglob_nb = self.ibool[ispec, ipt]
                    if node_mask[iglob_nb] != iglob:
                        self.node_adj_ptr[iglob + 1] += 1
                        node_mask[iglob_nb] = iglob

        # 3. Populate node adjacency data
        self.node_adj_data = np.zeros((self.node_adj_ptr[-1]), dtype=int)
        node_mask[:] = -1
        for iglob in range(self.nglob):
            istart = node2elem_ptr[iglob]
            iend = node2elem_ptr[iglob + 1]
            ic = 0
            for idx in range(istart, iend):
                ispec = node2elem_data[idx]
                for ipt in range(self.ngll):
                    iglob_nb = self.ibool[ispec, ipt]
                    if node_mask[iglob_nb] != iglob:
                        idx_nb = self.node_adj_ptr[iglob] + ic
                        self.node_adj_data[idx_nb] = iglob_nb
                        node_mask[iglob_nb] = iglob
                        ic += 1

    def setup_petsc(self):
        # Create l2g map and local node adjacency
        self._create_l2g_map()
        self._create_node_adjacency()

        comm = MPI.COMM_WORLD

        # We will use an array of Python sets to automatically filter duplicate global edges!
        # owned_edges[0] corresponds to local row_id 0, etc.
        owned_edges = [set() for _ in range(self.nglob_owned)]

        # 1. ADD LOCAL CONNECTIONS TO THE SETS
        for iglob in range(self.nglob):
            if self.owner_ranks[iglob] == self.myrank:
                row_id = self.l2g_map[iglob] - self.rstart
                istart = self.node_adj_ptr[iglob]
                iend = self.node_adj_ptr[iglob + 1]
                
                # Add all local neighbors (converted to global IDs) into the set
                for idx in range(istart, iend):
                    global_neighbor = self.l2g_map[self.node_adj_data[idx]]
                    owned_edges[row_id].add(global_neighbor)

        # 2. PACK THE GHOST EDGE PAIRS: [target_global_id, neighbor_global_id, ...]
        send_buffers = []
        size_sd = np.zeros(self.num_neighbors, dtype=np.int32)
        
        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            npts = self.neighbor_points[ib]
            
            edges_for_this_neighbor = []
            for ipt in range(npts):
                iglob = self.send_points[ib, ipt]
                
                # If they own it, I must send my connections to them
                if self.owner_ranks[iglob] == nrank:
                    target_global = self.l2g_map[iglob]
                    istart = self.node_adj_ptr[iglob]
                    iend = self.node_adj_ptr[iglob + 1]
                    
                    for idx in range(istart, iend):
                        neighbor_global = self.l2g_map[self.node_adj_data[idx]]
                        edges_for_this_neighbor.append(target_global)
                        edges_for_this_neighbor.append(neighbor_global)
            
            send_array = np.array(edges_for_this_neighbor, dtype=np.int64)
            send_buffers.append(send_array)
            size_sd[ib] = len(send_array)

        # 3. EXCHANGE SIZES
        size_rv = np.zeros(self.num_neighbors, dtype=np.int32)
        req_sd = []
        req_rv = []
        
        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            # mpi4py syntax: pass the numpy array slice to enforce fast C-level buffer passing
            req_sd.append(comm.Isend(size_sd[ib:ib+1], dest=nrank, tag=11))
            req_rv.append(comm.Irecv(size_rv[ib:ib+1], source=nrank, tag=11))
            
        MPI.Request.Waitall(req_rv) 
        MPI.Request.Waitall(req_sd)
        
        # 4. EXCHANGE THE EDGE LISTS
        recv_buffers = []
        req_sd = []
        req_rv = []
        
        for ib in range(self.num_neighbors):
            nrank = self.neighbor_ranks[ib]
            
            # Allocate exact sized receive buffer
            recv_array = np.zeros(size_rv[ib], dtype=np.int64)
            recv_buffers.append(recv_array)
            
            req_sd.append(comm.Isend(send_buffers[ib], dest=nrank, tag=12))
            req_rv.append(comm.Irecv(recv_array, source=nrank, tag=12))
            
        MPI.Request.Waitall(req_rv) 
        MPI.Request.Waitall(req_sd)

        # 5. MERGE INCOMING GHOST EDGES INTO LOCAL SETS
        for ib in range(self.num_neighbors):
            recv_array = recv_buffers[ib]
            
            # Unpack the pairs [target, neighbor, target, neighbor...]
            for i in range(0, len(recv_array), 2):
                target_global = recv_array[i]
                neighbor_global = recv_array[i+1]
                
                # We own the target, so calculate its local row index
                row_id = target_global - self.rstart
                
                # The 'set' automatically ignores it if we already found this connection locally!
                owned_edges[row_id].add(neighbor_global)

        # 6. COUNT d_nnz AND o_nnz
        d_nnz = np.zeros(self.nglob_owned, dtype=np.int32)
        o_nnz = np.zeros(self.nglob_owned, dtype=np.int32)
        
        for row_id in range(self.nglob_owned):
            for global_neighbor in owned_edges[row_id]:
                # Is the connection inside my owned block?
                if self.rstart <= global_neighbor < self.rend:
                    d_nnz[row_id] += 1
                else:
                    o_nnz[row_id] += 1

        self.d_nnz = d_nnz
        self.o_nnz = o_nnz

        # DEBUG PRINT
        mpi_print(f"d_nnz = {self.d_nnz}, o_nnz = {self.o_nnz}")

def main():

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()
    if nprocs != 2:
        if myrank == 0:
            print("This test is designed for 2 ranks, but got ", nprocs)
        exit(1)

    mesh = Mesh(nspec=2, ngll=4)
    mesh.build()
    mesh.setup_petsc()

if __name__ == "__main__":
    main()