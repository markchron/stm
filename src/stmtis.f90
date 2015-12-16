! METIS interface API
      module stmapimts
      ! Graph (V,E) data structure consists of two arrays:
      ! vtxs, info. about the vertices
      ! total vertex number
      !integer                   :: nvtxs
      ! total connections number
      integer                   :: nedges
      ! iadj(nvtxs+1), the index into adjncy that is the beginning of
      ! the adjacency list of vertices
      !integer, pointer          :: iadj
      ! adjncy, the adjacency lists of the vertices
      !integer, pointer          :: adjncy
      ! number of balancing constraints. >=1
      ! If 'ncon' is the number of weights associated with each vertex,
      ! vwgt(n*ncon), stores 'ncon' consecutive entries for each vertex.
      integer                   :: ncon
      ! metis options
      integer, dimension(0:39)  :: opts
      ! vwgt(NULL), weight of the vertex
      integer, pointer          :: vwgt
      ! vsize(NULL), the amount of data that needs to be sent due to the
      ! ith vertex.  - communication cost
      ! The amount of data that needs to be sent (communication cost,
      ! vsize) is different from the weight of the vertex (computational 
      ! cost, vwgt).
      integer, pointer          :: vsize
      ! adjwgt(NULL), weight of the edges
      ! contains 2*m entries (m=edges no.) The weight of edge geocadj(j)
      ! is stored at location adjwgt(j)
      !integer, pointer          :: adjwgt
      real(8), pointer          :: tpwgts
      real(8), pointer          :: ubvec
      integer                   :: objval

      ! variables pass by main program
      integer                   :: metis_rb 
      contains
      !
      subroutine metis_init
      metis_rb = 0
      ncon = 1
      nullify(vwgt, vsize)
      !nullify(adjwgt)
      nullify(tpwgts, ubvec)
      ! set option before calling METIS function
      call METIS_SetDefaultOptions(opts)
      end subroutine metis_init
! PURPOSE:
! partitioning the graph
! ndomain : domains number
! b_dist(nvtxs) 
      subroutine metis_api(nvtxs, nadys, iadj, adjncy, adjwgt, ndomain, b_dist)
      use stmheader, only : STDD
      integer, intent(in)                               :: nvtxs
      integer, intent(in)                               :: nadys
      integer, dimension(nvtxs+1), intent(in)           :: iadj
      integer, dimension(nadys),intent(in)              :: adjncy
      real(STDD), dimension(nadys), intent(in)          :: adjwgt

      integer, intent(in)                               :: ndomain
      integer, dimension(nvtxs), intent(out)            :: b_dist

      call metis_init
      ! iadj-adjncy, the adjacency structure (CSR) of the graph
      ! multilevel recursive bisection 
      if(metis_rb .eq. 1) then
        call set_metis_options_rb
        call METIS_PartGraphRecursive(nvtxs, ncon,iadj, adjncy, vwgt, vsize,adjwgt,ndomain, tpwgts, ubvec, opts, objval, b_dist)
      else
      ! multilevel k-way partition
        call set_metis_options_kway
        call METIS_PartGraphKway(nvtxs, ncon, iadj, adjncy, vwgt, vsize,adjwgt,ndomain, tpwgts, ubvec, opts, objval, b_dist)
      endif
      end subroutine metis_api
! METIS options
      ! options : recursive bisection
      subroutine set_metis_options_rb
      opts(0) = 0   ! 0, multilevel recursive bisectioing
      ! 2, METIS_OPTION_CTYPE, the matching schme to be used during
      ! coarsening. 
      opts(2) = 0
      ! 3, METIS_OPTION_IPTYPE, algorithms during initial partition
      opts(3) = 0
      ! 4, METIS_OPTION_RTYPE, algorithm for refinement
      opts(4) = 0
      ! 5, METIS_OPTION_DBGLVL, progress/debugging information. BITWISE
      opts(5) = 1 
      !opts(5) = 0
      ! 6, METIS_OPTION_NITER, the number of iterations for the
      opts(6) = 10
      ! 7, METIS_OPTION_NCUTS, the number of different partitiongs that
      opts(7) = 1
      ! 8, METIS_OPTION_SEED, seed for the random number generator
      opts(8) = 0
      ! 9, METIS_OPTION_NO2HOP, 2-hop matching
      opts(9) = 0
      ! 16, METIS_OPTION_UFACTOR, the maximum allowed load imbalance
      opts(16) = 1
      ! 17, METIS_OPTION_NUMBERING
      opts(17) = 1
      
      opts(1) = 0   ! 0, Edge-cut minimization
                    ! 1, total communication volume minimization
                    ! 2, METIS_OBJTYPE_NODE
      opts(10) = 0
      opts(11) = 1
      opts(12) = 0
      opts(13) = 0
      opts(14) = 0
      opts(15) = 1
      end subroutine set_metis_options_rb
      
      subroutine set_metis_options_kway
      ! 0, METIS_OPTION_PTYPE 
      opts(0) = 1       ! 0, (default) Multilevel recursive bisectioning
                        ! 1, Multilevel k-way partitioning
! multilevel k-way partitioning minimize the results subdomain
! connectiveity graph, enforce contiguous partitions, minimize
! alternative objectives and as such it may be preferable than
! multilevel recursive bisection.
      ! 1, METIS_OPTION_OBJTYPE, minimum objective function
      opts(1) = 0       ! 0, (default) Edge-cut minimization
                        ! 1, Total communication volume minimization
                        ! 2, METIS_OBJTYPE_NODE
      ! 2, METIS_OPTION_CTYPE, the matching schme to be used during
      ! coarsening. 
      opts(2) = 1       ! 0, (default) random matching
                        ! 1, sorted heavy-edge matching
      ! 3, METIS_OPTION_IPTYPE, algorithms during initial partition
      opts(3) = 2       ! 0, Grows a bisection using a greedy strategy
                        ! 1, computes a bisection at random followed by
                        ! refinements
                        ! 2, derives a separator from an edge cut
                        ! 3, grows a bisection using a greedy node-based
                        ! strategy
      ! 4, METIS_OPTION_RTYPE, algorithm for refinement
      opts(4) = 0       ! 0, FM-based cut refinement
                        ! 1, Greedy-based cut and volume refinement
                        ! 2, Two side node FM refinement
                        ! 3, One-side node FM refinement
      ! 5, METIS_OPTION_DBGLVL, progress/debugging information. BITWISE
      opts(5) = 0       ! 0, (default) 
      !opts(5) = 1 ! prints various diagnostic|debug msg
      ! 6, METIS_OPTION_NITER, the number of iterations for the
      ! refinement algorithm at each stage of the uncoarsening process
      opts(6) = 10      ! 10, default
      ! 7, METIS_OPTION_NCUTS, the number of different partitiongs that
      ! it will compute. 
      opts(7) = 1       !(default) The final partitioning
                        ! is the one that achieves the best edgecut or
                        ! communication volume.
      ! 8, METIS_OPTION_SEED, seed for the random number generator
      opts(8) = 0
      ! 9, METIS_OPTION_NO2HOP, 2-hop matching, is very effective for
      ! graphs with power-law degree distributions
      opts(9) = 0       ! 0, (default) perform a 2-hop matching
                        ! 1, not perform 
      ! 10, METIS_OPTION_MINCONN, the partitioning routines should try
      ! to minimize the maximum degree of the subdomain graph, i.e., the
      ! graph in which each partition is a node, and edges connect
      ! subdomains with a shared interface
      opts(10) = 1      ! 0, (default) does not explicitly minimize the maximum
                        ! connectivity
                        ! 1, explicitly minimize the maximum
                        ! connectivity
      ! 11, METIS_OPTION_CONTIG, produce partitions that are contiguous.
      ! If the input graph is not connected, this option is ignored
      opts(11) = 1      ! 0, not force to be contiguous partitions
                        ! 1, forces contiguous partitions
      ! 16, METIS_OPTION_UFACTOR, the maximum allowed load imbalance
      ! among the partitions. A value of x indicates that the allowed
      ! load imbalance is (1+x)/1000. 
      ! load imbalance of jth constraint = max_i(w(j,i)/t(j,i))
      ! w(j,i), the fraction of the overall weight of the jth constraint
      ! that is assigned to the ith partition 
      ! t(j,i), the desired target weigth of the jth constraint for the
      ! ith partition (specified by - tpwgts)
      opts(16) = 30     ! 1 for RB, 30 for KWAY
      ! 17, METIS_OPTION_NUMBERING, indicate which numbering scheme for
      ! adjacency structure of a graph or the element node structure of
      ! a mesh
      opts(17) = 1      ! 0, C-style numbering
                        ! 1, Fortran-style numbering 
!-- kway valid options
! --
      ! 12, METIS_OPTION_COMPRESS, graph is compressed by combining
      ! together vertices that have identical adjacency list
      opts(12) = 0      ! 0, not try to compress the graph
                        ! 1, compress the graph
      ! 13, METIS_OPTION_CCORDER,  if the connected components of the
      ! graph should first be identified and ordered separately
      opts(13) = 0      ! 0, not identify the connected components
                        ! 1, 
      ! 14, METIS_OPTION_PFACTOR, the miminum degree of the vertices
      ! that will be ordered last.
      opts(14) = 0      ! default, no vertices are removed
      ! 15, METIS_OPTION_NSEPS, the number of different separators that it
      ! will compute at each level of nested dissection. The final
      ! separator that is used is the smallest one
      opts(15) = 1      ! default
      end subroutine set_metis_options_kway
      end module stmapimts
