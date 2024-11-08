/***************************************************************************
 *                                                                         *
 *  DynELA Project                                                         *
 *                                                                         *
 *  (c) Copyright 1997-2004                                                *
 *                                                                         *
 *      Equipe C.M.A.O                                                     *
 *      Laboratoire Genie de production                                    *
 *      Ecole Nationale d'Ingenieurs de Tarbes                             *
 *      BP 1629 - 65016 TARBES cedex                                       *
 *                                                                         *
 *                                                                         *
 *  Main Author: Olivier PANTALE                                           *
 *                                                                         *
 **************************************************************************/

// begin date : 
// revision date : 
#include <femLibrary.h>
#include <InputDyn.h>
#include <lsdynaReader.h>
#include <femLibrary.h>

#include <Select.h>

#include <omp.h>


#include <VtkInterface.h>

#include <iostream>
#include <vector>
#include <set>


#include <iostream>

#include <vector>
#include <set>
#include <utility>
#include <utility>

#include <map>

/** Include the mmg2d library hader file */
// if the header file is in the "include" directory
// #include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2d.h"

String parsedFileName;

/* FROM INTERNAL
 * 
 * std::vector <std::pair<int,int>> getAllEdges(Structure *st){
  
      // Define the mesh as a vector of triangles, each represented by 3 vertex indices
    std::vector<std::vector<int>> triangles = {
        {1, 2, 3},
        {3, 2, 4},
        {4, 2, 5},
        {5, 2, 1}
    };

    // Map to store edges and their counts
    std::map<std::pair<int, int>, int> edgeCount;
    std::vector<std::pair<int, int>>  extEdges; //If we wnat only ext edges

    // Iterate over each triangle and its edges
    for (const auto& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            int v1 = triangle[i];
            int v2 = triangle[(i + 1) % 3]; // Next vertex in the triangle
            std::pair<int, int> edge = createEdge(v1, v2);
            edgeCount[edge]++;
        }
    }

    // Print the exterior edges
    std::cout << "Exterior edges:" << std::endl;
    for (const auto& edgeEntry : edgeCount) {
        if (edgeEntry.second == 1) { // Exterior edge appears only once
            std::cout << "Edge between vertices " << edgeEntry.first.first
                      << " and " << edgeEntry.first.second << std::endl;
           extEdges.push_back(edgeEntry.first); 
        }
    }  
  
  
}

----*/

// A helper function to create a sorted pair for an edge
std::pair<int, int> createEdge(int a, int b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::vector <std::pair<int,int>> getAllEdges(Structure *st, int el_ind_max){
 
 
    // Map to store edges and their counts
    std::map<std::pair<int, int>, int> edgeCount;
    std::vector<std::pair<int, int>>  extEdges; //If we wnat only ext edges


    int emax = 0;   
    int ecount = st->getElementsNumber();
    cout << "checking on "<< ecount <<" elements"<<endl;
    for (int e=0;e<el_ind_max;e++){
      Element *pel = st->getElement(e);
      
      for (int k = 0; k < pel->getNumberOfSideFaces(); k++){
        //int toHave = pel->getNumberOfNodesOnSideFace(k);
        //THIS IS LIKE IN SIDE.C
        Node* n1 = pel->getNodeOnSideFace(k, 0); //ASUMMIN
        Node* n2 = pel->getNodeOnSideFace(k, 1); //ASUMMIN

        //cout << "Elem "<< e <<" N1 "<<n1->Id<<", "<<n2->Id<<endl;
        std::pair<int, int> edge = createEdge(n1->Id, n2->Id); //COMPARE WITH number
        if (edge.second > emax) emax = edge.second ;
        edgeCount[edge]++; //THIS INSERT ELEMENT 
        
        //cout << "Edge count size "<<edgeCount.size()<<endl;
        //Node * 
        //Node *Element::getNodeOnEdge(short edge, short node)
        
        }
      
    }//Element
    cout << "DONE, max edge node index"<<emax<<endl;

    // Print the exterior edges
    std::cout << "Exterior edges:" << std::endl;
    for (const auto& edgeEntry : edgeCount) {
        if (edgeEntry.second == 1) { // Exterior edge appears only once
            //std::cout << "Edge between vertices " << edgeEntry.first.first
              //        << " and " << edgeEntry.first.second << std::endl;
           extEdges.push_back(edgeEntry.first); 
        }
    }  
  
  cout << "Edge sizes All:" <<edgeCount.size()<<", External: "<<extEdges.size()<<endl;
  return extEdges;
}

int main() {

    //Structure model;
    
  Domain *model = new Domain();
  Global_Structure = new Structure();


  /*
  model->createNode(1, 0.,0.,0.);
  model->createNode(2, .1,0.,0.);
  model->createNode(3, 0.,.1,0.);
  model->createNode(4, .1,.1,0.);


  model->createNode(5, 0.5,1.,0.);
  model->createNode(6, 1.5,1.,0.);
  model->createNode(7, 0.5,2.,0.);
  model->createNode(8, 1.5,2.,0.);
*/


  
  Real elem_x, elem_y;
  Element* el = new ElQua4nAx(1);

/*
  Real nbNodes=1;
  Real i;
  Real j;
  for (j=0;j<=nbElementsHauteur;j+=1) 
    for (i=0;i<=nbElementsLargeur;i+=1) {
      struct.createNode(nbNodes,i*dxLargeur,j*dxHauteur,0);
      cylinderNds.add(nbNodes);
      nbNodes++;
  };
  nbNodes--;
  */

  Element* el2 = new ElQua4nAx(2);
  Indice *ind = new Indice[4];
  ind[0]=1;ind[1]=2;ind[2]=4;ind[3]=3;

  Indice *ind2 = new Indice[4];
  ind2[0]=5;ind2[1]=6;ind2[2]=8;ind2[3]=7;  

  //model->add(el);
  //model->add(el2);
  
  //model->createElement(el,ind);
  cout << "Elem size "<<model->elements.size()<<endl;
  //model->createElement(el2,ind2);
  
  Global_Structure->setDefaultElement(el);
  //Global_Structure->createElement(1, 1, 2, 4, 3);
  

  // ElementSet allES("ES_All");
  // model.add(&allES,1,1);

  ElementSet allES();
  
  //lsdynaReader("sphere-plate.k");
  // model.add(&allES,1,1);

  //ASSIGN DOMAIN TO STRUCTURE
 //BEFORE NODESET CREATION!
 //Global_Structure->domains.add(model);
 //Global_Structure->setDomain(model);
 Global_Structure->domains(0)=model;
  NodeSet topNS;

  //topNS.add(3);
  //topNS.add(4);
  


  BoxMesher bm;
  bm.link(Global_Structure);
  bm.setElement(el);
  
  
  
  //bm.elementType = elementType
  
  bm.rectangle (1.0,1.0,10,10);

  cout << "Elements created: "<<Global_Structure->getElementsNumber()<<endl;
  cout << "Nodes created: "<<Global_Structure->getNodesNumber()<<endl;
  cout << "ADDING NODE ELEMENTS (NEEDED FOR CONTACT)"<<endl;

  Real i;
  Real j;
    
  j = 10;
  
  for (int i=0;i<=10;i++) {
    cout << "i j "<<i<<", "<<j<<endl;
    int nind = 11*j+i;
    int eind = 9*j+i;
    cout << "node ind , el ind "<<nind<<","<<eind<<endl;
    //if (i>0) Global_Structure->getNode(i)->elements<<Global_Structure->getElement(j*10+i);
    if (i<10) Global_Structure->getNode(nind)->elements<<Global_Structure->getElement(eind);
    if (i>0)  Global_Structure->getNode(nind)->elements<<Global_Structure->getElement(eind-1);
    
    cout << "Node "<< 10*j+i << "element size"<<Global_Structure->getNode(nind)->elements.size()<<endl;
  }
  cout << "FIRST ELEMENT INDEX AFTER FIRST BOX "<<Global_Structure->getElementsNumber()+1<<endl; 
          
  cout << "Done "<<endl;

  int nbNodes=Global_Structure->getNodesNumber()+1;
  int nel = Global_Structure->getElementsNumber();

  int nn = nbNodes; //BEFORE CREATION
  
  cout << "Node count before top mesh "<<nbNodes<<endl;
  for (j=0;j<=1;j+=1) 
    for (i=0;i<=10;i+=1) {
      //cout <<"x "<<i*0.1+0.5<<endl;
      //cout << "y "<<j*0.1+1.0<<endl;
      Global_Structure->createNode(nbNodes,i*0.045,j*0.10+1.00001,0);
      //cylinderNds.add(nbNodes);
      nbNodes++;
  };


  
  //nbNodes--;
  nn--;
  
  cout << "Nodes created: "<<Global_Structure->getElementsNumber()<<endl;
  
  cout << "Connectivity"<<endl;
  int nx = 10.0;
  Global_Structure->setDefaultElement(el);
  Indice nbElements = Global_Structure->getElementsNumber() + 1;
  Indice x1;
  Indice x2;
  Indice x3;
  Indice x4;

  for (Indice j = 0; j < 1; j++){
    for (Indice i = 0; i < 10; i++)
    {
      Global_Structure->setDefaultElement(el);  
      cout << "i "<<i<<endl;
      x1 = nn+(i + (j * (nx + 1)) + 1);
      x2 = nn+(i + (j * (nx + 1)) + 2);
      x3 = nn+(i + ((j + 1) * (nx + 1)) + 2);
      x4 = nn+(i + ((j + 1) * (nx + 1)) + 1);
      cout << "x1--x4 "<<x1 <<","<<x2 <<","<<x3 <<","<<x4 <<endl;
      Global_Structure->createElement(nbElements, x1, x2, x3, x4);
      nbElements++;
    }
  }      

  j = 10;
  for (int i=0;i<=10;i++) {
    
    int nind = nn + i;
    int eind = nel+ i;

    cout << "CREATING NODE ELEMENTS FOR MASTER "<<endl;
    cout << "node ind , el ind "<<nind<<","<<eind<<endl;
    if (i<10) Global_Structure->getNode(nind)->elements<<Global_Structure->getElement(eind);
    //if (i>0)  Global_Structure->getNode(nind)->elements<<Global_Structure->getElement(eind-1);
  }
  
    
  //bm2.rectangle (1.0,0.1,10,1);  
  
  cout << "Elements created: "<<Global_Structure->getElementsNumber()<<endl;
  cout << "Nodes created: "<<Global_Structure->getNodesNumber()<<endl;
  //cout << "Side elsize size: " << elside.sides(0)->size()<<endl;
  
  // omp_set_num_threads(1);

   nbNodes = Global_Structure->getNodesNumber();
  cout << "POINTS " << nbNodes << " float\n";

  for (long i = 0; i < nbNodes; i++)
    cout << Global_Structure->getNode(i)->coords(0)
            << " " << Global_Structure->getNode(i)->coords(1) << " "
            <<Global_Structure->getNode(i)->coords(2) << "\n";


   cout << "OVERALL ELEMENT NUMBER" << Global_Structure->getElementsNumber()<<endl;

  Interface int_body(1);
  
  Side masterside;
  //SideFace2D *sf = new SideFace2D;
  //masterside.addSideFace(sf);
  
  Side slaveside;
  
  //sf->addNode(model->getNodeByNumber(3));  
  //sf->addNode(model->getNodeByNumber(4));
  //elside.addSideFace(sf);

  NodeSet masterNS;
  NodeSet slaveNS;
  //start, end, inc (def 1)
  slaveNS.add(111,121);


  masterNS.add(122,132); //(TOP)

  masterside.addNodeSet(&masterNS);
  cout << "MASTER NODE SIZE "<<masterside.nodes.size()<<endl;

  slaveside.addNodeSet(&slaveNS);

  cout << "Init Side SLAVE "<<endl;


  CoulombLaw* cl = new CoulombLaw();
  cl->setFriction(0.0);
  int_body.contactLaw = cl;
  if (int_body.setMaster(&masterside) !=Success) {cout << "set master FAILS"<<endl;}
  int_body.setSlave(&slaveside);

  int_body.Init();
  Global_Structure->addInterface(&int_body);
  //////////////// CONTACT
  //Global_Structure->addInterface(&int_body);

  cout << "---------------------------\n Node "<< 1 << "element size"<<Global_Structure->getNode(1)->elements.size()<<endl;  
 
  Material steel;

  Real A=300e6;
  Real B=100e6;
  Real n=1;
  Real young=117e9;
  Real poisson= .35;
  Real density= 2700.0;

  IsoHardElastoplastic hard;
  hard.setYieldStress(A);
  hard.setHardParameter(B);
  hard.setHardExponent(n);

  steel.setYoung(young);
  steel.setPoisson(poisson);
  steel.setDensity(density);
  steel.setColor(1,0,1);
  steel.setHardening(&hard);
//steel.

steel.setHeatCoefficient(3.7900000E+02);
steel.setDilatation(1.2000000E-05);
steel.setInitTemperature(3.0000000E+02);
steel.setConductivity(4.6000000E+01);

  cout << "DOMAIN SIZE"<<Global_Structure->domains.size()<<endl;
  ElementSet all;
  for (int e=0;e<Global_Structure->getElementsNumber();e++){
    all.add(Global_Structure->domains(0)->elements(e));
    }
  
  Global_Structure->attachMaterialToElements(&steel,&all);
  
  //el->attachMaterial(&steel); //OTHERWISE CRASHES
  //el2->attachMaterial(&steel); //OTHERWISE CRASHES
  // steel.setHardeningLaw(hardLaw);
  // steel.youngModulus = young;
  // steel.poissonRatio = poisson;
  // steel.density = density;

  // model.add(&steel, &allES);




    NodeSet base;
    base.add(1,11);

    BoundaryRestrain baseDisp;
    baseDisp.set(0, 1, 0);
    Global_Structure->attachConstantBC(&baseDisp,&base);

    
    NodeSet axis;
    axis.add(1,143,11);
    BoundaryRestrain axisDisp;
    axisDisp.set(1, 0, 0);
    Global_Structure->attachConstantBC(&axisDisp,&axis);
    
    
    
    NodeSet top;
    top.add(132,143);


    BoundarySpeed topSpeed;
    topSpeed.set(0.0000000E+00, -5.0, 0.0000000E+00);
    Global_Structure->attachConstantBC(&topSpeed,&top);

  


  Global_Structure->name = "test";
  
  //Global_Structure->domains.add(model);
  
  //cout << "Domain size "<< Global_Structure->domains.size();
  //Global_Structure->setDomain(model);
  
  //ExplicitSolverCH *solver = new ExplicitSolverCH();
  

  ExplicitSolver *solver = new ExplicitSolver();
  

  solver->setTimeStepMethod("Courant");
   
  Global_Structure->addSolver(solver);
  Global_Structure->logFile = new LogFile("log.txt");
  
  //model->solvers.add(solver);

  //model->initSolve();
  //model->currentSolver=solver;
  //model->solve();


 omp_set_num_threads(4);
 cout << "THREADS "<<omp_get_max_threads<<endl;


 
  //Real stopTime=80.0e-6;
  Real stopTime=5.0e-2;
  Real saveTime=stopTime/40.0;
  solver->setTimes(0.0,stopTime);
  
  
  Global_Structure->resultFile=new io_Data;
  Global_Structure->resultFile->name = "test";
  Global_Structure->resultFile->pstructure = Global_Structure;
  Global_Structure->resultFile->setMode(Write);
  



  Global_Structure->setSaveTimes(0,stopTime,saveTime);

  //Global_Structure->setResultFile("results.bin");
    
  Global_Structure->initSolve();

  Global_Structure->vtk = new VtkInterface ("out.vtk");
  Global_Structure->vtk->pstructure = Global_Structure;
  
  //cout << "STRUCT ELEMENTS "<<Global_Structure->elements.size()<<endl;
  Global_Structure->solve();

  
  cout << "REMESHING ..."<<endl;
  ////// MMG THINGS
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;

  mmgMesh = NULL;
  mmgSol  = NULL;
  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);



  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** read the mesh in a mesh file */
  //MMG2D_loadMesh(mmgMesh,filename);
  
  
  std::vector <std::pair<int,int>> ext_edges = getAllEdges(Global_Structure, 100);
  
  //DEFINED IN :   mmg2d/API_functions_2d.c
  //MMG2D_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int nt, MMG5_int nquad, MMG5_int na);
  //int np    = Global_Structure->getNodesNumber();
  int np = 121;
  int quad, nt;
  int na    = ext_edges.size(); //EDGES
  //na: Number of edges
 
  bool split_quads = true;
  
  if (!split_quads){
  nquad = Global_Structure->getElementsNumber();
  nt    = 0; //TRIS    
    
  } else {
    nquad = 0;
    nt = 2 * Global_Structure->getElementsNumber();
    
  }
 
 //https://github.com/tan2/DynEarthSol/blob/master/remeshing.cxx
 
  
  //In API_functions
  //int MMG2D_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int nt, MMG5_int nquad, MMG5_int na) {
  if (MMG2D_Set_meshSize(mmgMesh, np,  nt,  nquad, na)==0)
    cout << "ERROR ALLOCATING MESH"<<endl;
  else 
    cout << "MESH CREATED OK"<<endl;
  cout << "Number of points: "<< mmgMesh->na << endl;

  



  int *edges = new int[2 * na]; 

  //int MMG2D_Set_vertex(MMG5_pMesh mesh, double c0, double c1, MMG5_int ref, MMG5_int pos)
  for (int n=0;n<121;n++){
    if (!MMG2D_Set_vertex(mmgMesh, Global_Structure->getNode(n)->coords(0), Global_Structure->getNode(n)->coords(1), NULL, n))
      cout << "ERROR ALLOCATING NODE "<<n<<endl;
  }
  cout << "Vertices allocated"<<endl;
  //int *ref = new int[na];
  ////// MMG2D_Set_edges IS CRASHING
  //int MMG2D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int *refs)
  //int res = MMG2D_Set_edges(mmgMesh, edges, nullptr);
  for (int e=0;e<na;e++)
    if (MMG2D_Set_edge(mmgMesh, ext_edges[e].first, ext_edges[e].second, NULL, e) !=1)
      cout << "ERROR CREATING EDGE "<<endl;

  cout << "EDGES ALLOCATED "<<endl;

  if (!split_quads){
  //int MMG2D_Set_quadrilateral(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int v3, MMG5_int ref, MMG5_int pos)
  //int  MMG2D_Set_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quadra, MMG5_int *refs) {
  for (int e=0;e<100;e++){
    MMG2D_Set_quadrilateral(mmgMesh, Global_Structure->getElement(e)->nodes(0)->Id
                                   , Global_Structure->getElement(e)->nodes(1)->Id
                                   , Global_Structure->getElement(e)->nodes(2)->Id
                                   , Global_Structure->getElement(e)->nodes(3)->Id
                                   , NULL, e);
  
  
  }
  
  } else { //split quads = true
    //int  MMG2D_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs)
    //int MMG2D_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int ref, MMG5_int pos)
    for (int e=0;e<200;e+=2){
    MMG2D_Set_triangles(mmgMesh,  Global_Structure->getElement(e)->nodes(0)->Id,
                                  Global_Structure->getElement(e)->nodes(1)->Id,
                                  Global_Structure->getElement(e)->nodes(2)->Id,
                                  NULL, e);
  
    MMG2D_Set_triangles(mmgMesh,  Global_Structure->getElement(e)->nodes(1)->Id,
                                  Global_Structure->getElement(e)->nodes(2)->Id,
                                  Global_Structure->getElement(e)->nodes(3)->Id,
                                  NULL, e+1);
    }
  }
  /** Set parameters : for example set the maximal edge size to 0.1 */
  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,0.1);

  /** Higher verbosity level */
  MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,5);


  /////// Generate the mesh ///////
  int ier = MMG2D_mmg2dmesh(mmgMesh,mmgSol);
  
  cout << "New mesh npoints "<<mmgMesh->np<<endl;

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DMESH: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DMESH\n");

  /*save result*/
  //if ( MMG2D_saveMesh(mmgMesh,outname) != 1 )
  //  exit(EXIT_FAILURE);

  /*save metric*/
  //if ( MMG2D_saveSol(mmgMesh,mmgSol,outname) != 1 )
  //  exit(EXIT_FAILURE);


  
  
  Global_Structure->delAllData();
  
  
  
  
  Domain *dom = new Domain();
  Global_Structure->setDomain( dom);
  cout << "CURRENT DOMEL  SIZE "<<Global_Structure->getCurrentDomain()->elements.size()<<endl;
  
  //Global_Structure->createNode(0,0,0,0);
  
  cout << "struct nodecount "<<Global_Structure->getNodesNumber()<<endl;
  
  
  cout << "Node count" << mmgMesh->np<<endl;
  cout <<"Mesh node 0 "<<mmgMesh->point[0].c[0]<<endl;
  
  for (int n=0;n<mmgMesh->np;n++)
    Global_Structure->createNode(n, mmgMesh->point[n].c[0], mmgMesh->point[n].c[1], 0);


  Element* el3 = new ElTri3n2D();
  Global_Structure->setDefaultElement(el3);
  
  for (int tri=0;tri<mmgMesh->nt;tri++){
    //cout << "\ntria "<<tri<<endl;
    //cout << mmgMesh->tria[tri].v[0] <<", "<<
    //mmgMesh->tria[tri].v[1] <<", "<<
    //mmgMesh->tria[tri].v[2] <<", "<<"NP"<<mmgMesh->np<<endl;
    bool error = false;
    for (int i=0;i<3;i++){ 
//      cout << "i "<<i<<endl;
      if (mmgMesh->tria[tri].v[i] >= Global_Structure->getNodesNumber()){
        cout << "ERROR on INDEX "<<endl;
        error = true;
      }
    }
    if (!error)
    Global_Structure->createElement(tri,mmgMesh->tria[tri].v[0] ,
                                        mmgMesh->tria[tri].v[1] , 
                                        mmgMesh->tria[tri].v[2] );
                                        
  }
  
  VtkInterface out;
  out.openFile("test.vtk");
  //out.dataWrite();
  out.pstructure = Global_Structure;  
  out.headerWrite();
  out.nodesWrite();

  // Write the nodes
  out.elementsWrite();

  
  out.close();

/*
  int ge = 0; //global elem
  Global_Structure->setDefaultElement(el);
  cout << "\nReallocating mesh" <<endl;
  for (int q=0;q<mmgMesh->nquad;q++){
    cout << "quad "<<q<<endl;
    
      Global_Structure->createElement(ge, mmgMesh->quadra[q].v[0], 
                                                  mmgMesh->quadra[q].v[1],
                                                  mmgMesh->quadra[q].v[2],
                                                  mmgMesh->quadra[q].v[3]);    
    ge++;
  }
*/  
  //cout << "Struct El count "<<Global_Structure->getElementsNumber()<<endl;
     
        
  //inout_s.c
  
  /** Set parameters : for example set the maximal edge size to 0.1 */
  //MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,0.1);

  /** Higher verbosity level */
  //MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,5);


  /** Generate the mesh */
  //ier = MMG2D_mmg2dmesh(mmgMesh,mmgSol);


/*
  // Set parameters : for example set the maximal edge size to 0.1 
  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,0.1);

  // Higher verbosity level
  MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,5);


  // Generate the mesh 
  ier = MMG2D_mmg2dmesh(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DMESH: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DMESH\n");
*/

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);
  
  
  //vtk.write();
    
  return 0;
  
  
  
}
