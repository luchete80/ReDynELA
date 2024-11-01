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
  Element* el = new ElQua4nAx(1);
/*
  Element* el2 = new ElQua4nAx(2);
  Indice *ind = new Indice[4];
  ind[0]=1;ind[1]=2;ind[2]=4;ind[3]=3;

  Indice *ind2 = new Indice[4];
  ind2[0]=5;ind2[1]=6;ind2[2]=8;ind2[3]=7;  
*/
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
  
  Interface int_body(1);
  
  Side elside;
  SideFace2D *sf = new SideFace2D;
  
  //sf->addNode(model->getNodeByNumber(3));  
  //sf->addNode(model->getNodeByNumber(4));
  //elside.addSideFace(sf);
	
  elside.addNodeSet(&topNS);
  cout << "Init Side "<<endl;
  elside.Init();


  BoxMesher bm;
  bm.link(Global_Structure);
  bm.setElement(el);
  //bm.elementType = elementType
  
  bm.rectangle (1.0,1.0,10,10);
  
  
  //cout << "Side elsize size: " << elside.sides(0)->size()<<endl;
  
  // omp_set_num_threads(1);
	
  //NodeSet bottomNS;   
  //bottomNS.add(1,2,1);
  // NodeSet bottomNSy("NS_Bottomy"); model.add(&bottomNS, 2);
  // NodeSet bottomNSx("NS_Bottomx"); model.add(&bottomNS, 3);
  // NodeSet bottomNSz("NS_Bottomz"); model.add(&bottomNS, 4);
 
  // ElasticLaw *hardLaw = new ElasticLaw;
 
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


  // BoundaryRestrain bottomBC ("BC_bottom");
  // bottomBC.setValue(1, 1, 1);
  // model.attachConstantBC(&bottomBC, &bottomNS);

  // BoundaryRestrain bottomBCx ("BC_bottomxz");
  // bottomBC.setValue(1, 0, 1);
  // model.attachConstantBC(&bottomBCx, &bottomNS);
  
  // BoundaryRestrain bottomBCy ("BC_bottomyz");
  // bottomBC.setValue(0, 1, 1);
  // model.attachConstantBC(&bottomBCy, &bottomNS);

  // BoundaryRestrain bottomBCz ("BC_bottomz");
  // bottomBC.setValue(0, 0, 1);
  // model.attachConstantBC(&bottomBCx, &bottomNS);  

  // BoundarySpeed speedBC ("BC_speed");
  // speedBC.setValue(0, 0, -speed);
  // model.attachConstantBC(&speedBC, &topNS);

  // Explicit solver("Solver");
  // solver.setTimes(0, stopTime);
  // model.add(&solver);
  // model.setSaveTimes(0, stopTime, stopTime / nbreSaves);


  // HistoryFile vonMisesHist("vonMisesHistory");

  // vonMisesHist.setFileName("vonMises.plot");
  // vonMisesHist.add(&allES, 0, Field::vonMises);

  // model.add(&vonMisesHist);


  // model.solve();
  
  
  
  //pdomain->add(el); //THIS CRASHES
  //pdomain->typeOfPreprocessing=preprocessingSingle;

  // initialisation des donnees aux noeuds
  //model->initSolve();
  
  Global_Structure->name = "test";
  
  //Global_Structure->domains.add(model);
  
  //cout << "Domain size "<< Global_Structure->domains.size();
  //Global_Structure->setDomain(model);
  
  ExplicitSolver *solver = new ExplicitSolver();

  //ExplicitSolver *solver = new ExplicitSolverCH();
  
  solver->setTimes(0.0,1.0);
  solver->setTimeStepMethod("Courant");
   
  Global_Structure->addSolver(solver);
  Global_Structure->logFile = new LogFile("log.txt");
  
  //model->solvers.add(solver);

  //model->initSolve();
  //model->currentSolver=solver;
  //model->solve();


 omp_set_num_threads(4);
 cout << "THREADS "<<omp_get_max_threads<<endl;
 
Real stopTime=80.0e-6;

  Global_Structure->resultFile=new io_Data;
  Global_Structure->resultFile->name = "test";
  Global_Structure->resultFile->pstructure = Global_Structure;
  Global_Structure->resultFile->setMode(Write);


  Global_Structure->setSaveTimes(0,stopTime,stopTime/1.0);

  //Global_Structure->setResultFile("results.bin");
    
  Global_Structure->initSolve();

  //cout << "STRUCT ELEMENTS "<<Global_Structure->elements.size()<<endl;
  Global_Structure->solve();
  
  return 0;
  
  
  
}
