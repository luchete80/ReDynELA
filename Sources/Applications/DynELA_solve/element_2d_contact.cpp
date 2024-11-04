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

String parsedFileName;

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

  

  Global_Structure->setDefaultElement(el);  
  int nbNodes=Global_Structure->getNodesNumber()+1;
  int nn = nbNodes;
  Real i;
  Real j;
  for (j=0;j<=1;j+=1) 
    for (i=0;i<=10;i+=1) {
      cout <<"x "<<i*0.1+0.5<<endl;
      cout << "y "<<j*0.1+1.0<<endl;
      Global_Structure->createNode(nbNodes,i*0.1,j*0.10+1.00001,0);
      //cylinderNds.add(nbNodes);
      nbNodes++;
  };
  //nbNodes--;

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
    for (Indice i = 0; i < 9; i++)
    {
      cout << "i "<<i<<endl;
      x1 = nn+(i + (j * (nx + 1)) + 1);
      x2 = nn+(i + (j * (nx + 1)) + 2);
      x3 = nn+(i + ((j + 1) * (nx + 1)) + 2);
      x4 = nn+(i + ((j + 1) * (nx + 1)) + 1);
      Global_Structure->createElement(nbElements, x1, x2, x3, x4);
      nbElements++;
    }
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




  Interface int_body(1);
  
  Side masterside;
  //SideFace2D *sf = new SideFace2D;
  //masterside.addSideFace(sf);
  
  Side slaveside;
  
  //sf->addNode(model->getNodeByNumber(3));  
  //sf->addNode(model->getNodeByNumber(4));
  //elside.addSideFace(sf);

  NodeSet masterNS;
  //start, end, inc (def 1)
  masterNS.add(111,121);

  NodeSet slaveNS;
  slaveNS.add(122,132);

  masterside.addNodeSet(&masterNS);
  slaveside.addNodeSet(&slaveNS);
  cout << "Init Side "<<endl;
  masterside.Init();
  slaveside.Init();
  CoulombLaw* cl = new CoulombLaw();
  cl->setFriction(0.2);
  int_body.contactLaw = cl;
  int_body.setMaster(&masterside);
  int_body.setSlave(&slaveside);
  //////////////// CONTACT
  //Global_Structure->addInterface(&int_body);
  
 
  Material steel;

  Real A=400e6;
  Real B=100e6;
  Real n=1;
  Real young=117e9;
  Real poisson= .35;
  Real density= 8930.0;

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
    axis.add(1,143,10);
    BoundaryRestrain axisDisp;
    axisDisp.set(1, 0, 0);
    Global_Structure->attachConstantBC(&axisDisp,&axis);
    
    
    
    NodeSet top;
    top.add(132,143);


    BoundarySpeed topSpeed;
    topSpeed.set(0.0000000E+00, -1.0, 0.0000000E+00);
    Global_Structure->attachConstantBC(&topSpeed,&top);

  


  Global_Structure->name = "test";
  
  //Global_Structure->domains.add(model);
  
  //cout << "Domain size "<< Global_Structure->domains.size();
  //Global_Structure->setDomain(model);
  
  //ExplicitSolverCH *solver = new ExplicitSolverCH();
  

  ExplicitSolver *solver = new ExplicitSolver();
  
  solver->setTimes(0.0,1.0e-1);
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
  Real stopTime=5.0e-5;

  Global_Structure->resultFile=new io_Data;
  Global_Structure->resultFile->name = "test";
  Global_Structure->resultFile->pstructure = Global_Structure;
  Global_Structure->resultFile->setMode(Write);
  



  Global_Structure->setSaveTimes(0,stopTime,stopTime);

  //Global_Structure->setResultFile("results.bin");
    
  Global_Structure->initSolve();

  Global_Structure->vtk = new VtkInterface ("out.vtk");
  Global_Structure->vtk->pstructure = Global_Structure;
  
  //cout << "STRUCT ELEMENTS "<<Global_Structure->elements.size()<<endl;
  Global_Structure->solve();


  
  //vtk.write();
    
  return 0;
  
  
  
}
