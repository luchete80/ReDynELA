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

    Structure model;

  model.createNode(1, 0.,0.,0.);
  model.createNode(2, .1,0.,0.);
  model.createNode(3, 0.,.1,0.);
  model.createNode(4, .1,.1,0.);
  model.createNode(5, 0.,0.,.1);
  model.createNode(6, .1,0.,.1);
  model.createNode(7, 0.,.1,.1);
  model.createNode(8, .1,.1,.1);

  Element* el = new ElHex8n3D;
  model.setDefaultElement(el);
  model.createElement(1, 1, 2, 4, 3, 5, 6, 8, 7);


  // ElementSet allES("ES_All");
  // model.add(&allES,1,1);

  ElementSet allES();
  
  lsdynaReader("sphere-plate.k");
  // model.add(&allES,1,1);

  NodeSet topNS;

  topNS.add(5);
  topNS.add(6);
  topNS.add(7);
  topNS.add(8);
  
	// omp_set_num_threads(1);
	
  NodeSet bottomNS;   
  bottomNS.add(1,4,1);
  // NodeSet bottomNSy("NS_Bottomy"); model.add(&bottomNS, 2);
  // NodeSet bottomNSx("NS_Bottomx"); model.add(&bottomNS, 3);
  // NodeSet bottomNSz("NS_Bottomz"); model.add(&bottomNS, 4);
 
  // ElasticLaw *hardLaw = new ElasticLaw;
 
  Material steel;
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

  
  return 0;
  
}