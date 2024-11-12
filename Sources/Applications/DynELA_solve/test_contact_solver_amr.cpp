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

#define MAX4(a,b,c,d)  (((MAX0(a,b)) > (MAX0(c,d))) ? (MAX0(a,b)) : (MAX0(c,d)))
#define MAX0(a,b)     (((a) > (b)) ? (a) : (b))

String parsedFileName;

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


  Real elem_x, elem_y;
  Element* el = new ElQua4nAx(1);

  Element* el2 = new ElQua4nAx(2);
  Indice *ind = new Indice[4];
  ind[0]=1;ind[1]=2;ind[2]=4;ind[3]=3;

  Indice *ind2 = new Indice[4];
  ind2[0]=5;ind2[1]=6;ind2[2]=8;ind2[3]=7;  

  
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
  

  ExplicitSolverRmsh *solver = new ExplicitSolverRmsh();
  

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
  
  Global_Structure->reMesh();
  
  
  
}
