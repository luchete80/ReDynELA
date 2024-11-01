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

String parsedFileName;

int main() {

    //Structure model;
    
  Domain *model = new Domain();
  Global_Structure = new Structure();
  
  model->createNode(1, 0.,0.,0.);
  model->createNode(2, .1,0.,0.);
  model->createNode(3, 0.,.1,0.);
  model->createNode(4, .1,.1,0.);



  Element* el = new ElQua4nAx;
  Indice *ind = new Indice[4];
  ind[0]=1;ind[1]=2;ind[2]=4;ind[3]=3;


  model->createElement(el,ind);

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

  topNS.add(3);
  topNS.add(4);
  
  Interface int_body(1);
	// omp_set_num_threads(1);
	
  //NodeSet bottomNS;   
  //bottomNS.add(1,2,1);
  // NodeSet bottomNSy("NS_Bottomy"); model.add(&bottomNS, 2);
  // NodeSet bottomNSx("NS_Bottomx"); model.add(&bottomNS, 3);
  // NodeSet bottomNSz("NS_Bottomz"); model.add(&bottomNS, 4);
 
  // ElasticLaw *hardLaw = new ElasticLaw;
 
  Material steel;
  
  el->attachMaterial(&steel); //OTHERWISE CRASHES
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
  
  Global_Structure->domains.add(model);
  cout << "Domain size "<< Global_Structure->domains.size();
  //Global_Structure->setDomain(model);
  
  ExplicitSolver *solver = new ExplicitSolver();
  Global_Structure->addSolver(solver);
  Global_Structure->logFile = new LogFile("log.txt");
  
  model->solvers.add(solver);
  
  model->initSolve();
  
  //Global_Structure->initSolve();
  
  //Global_Structure->solve();
  
  return 0;
  
  
  
}
