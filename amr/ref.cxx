//PORTION OF


    MMG2D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet,
                    &mmgSol, MMG5_ARG_end);

    // 2) Build mesh in MMG5 format
    // Manually set of the mesh 
    //  a) give the size of the mesh: vertices, triangles, quads(=0), edges
    if ( MMG2D_Set_meshSize(mmgMesh, old_nnode, old_nelem, 0, old_nseg) != 1 )
        exit(EXIT_FAILURE);
    //   b) give the vertex coordinates. References are NULL but can be an integer array for boundary flag etc.
    if( MMG2D_Set_vertices(mmgMesh, qcoord, NULL) != 1)
        exit(EXIT_FAILURE);
    //   c) give the connectivity. References are NULL but can be an integer array for boundary flag etc.
    for (int i = 0; i < old_nelem*NODES_PER_ELEM; ++i)
        ++qconn_from_1[i];
    if( MMG2D_Set_triangles(mmgMesh, qconn_from_1, NULL) != 1 )
        exit(EXIT_FAILURE);
    //   d) give the segments (i.e., boundary edges)
    for (int i = 0; i < old_nseg*NODES_PER_FACET; ++i)
        ++qsegment_from_1[i];
    for (int i = 0; i < old_nseg; ++i)        
        if( MMG2D_Set_edge(mmgMesh, qsegment_from_1[i*NODES_PER_FACET], 
                qsegment_from_1[i*NODES_PER_FACET+1], qsegflag[i], i+1) != 1)
            exit(EXIT_FAILURE);
    // if( MMG2D_Set_edges(mmgMesh, qsegment_from_1, qsegflag) != 1 )
    //     exit(EXIT_FAILURE);

    // 3) Build sol in MMG5 format
    //      Here a 'solution' is a nodal field that becomes
    //      the basis for metric tensor for isotropic and anisotropic
    //      mesh adaptation.
    //
    //   a) give info for the sol structure
    if( MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, old_nnode, MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
    //   b) give solutions values and positions
    compute_metric_field(var, old_connectivity, param.mesh.resolution, *var.ntmp, *var.tmp_result_sg);
    //      i) If sol array is available:
    if( MMG2D_Set_scalarSols(mmgSol, (*var.ntmp).data()) != 1 )
        exit(EXIT_FAILURE);
    //      ii) Otherwise, set a value node by node:
    // for (int i = 0; i < var.nnode; ++i) {
    //     if( MMG2D_Set_scalarSol(mmgSol, 0.5, i+1) != 1 )
    //         exit(EXIT_FAILURE);
    // }
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_optim, 0) != 1 )
    exit(EXIT_FAILURE);
