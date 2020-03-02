/***************************************************************************
 *                                                                         *
 *  DynELA Project                                                         *
 *                                                                         *
 *  (c) Copyright 1997-2006                                                *
 *                                                                         *
 *      Equipe C.M.A.O                                                     *
 *      Laboratoire Genie de production                                    *
 *      Ecole Nationale d'Ingenieurs de Tarbes                             *
 *      BP 1629 - 65016 TARBES cedex                                       *
 *                                                                         *
 *                                                                         *
 *  Main Author: Olivier PANTALE                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 **************************************************************************/

// begin date : 13/03/1997

/*
  Class ElQua4nAx definition
*/

#ifndef __ElQua4nAx_h__
#define __ElQua4nAx_h__

#define Name_ElQua4nAx "ElQua4nAx"

class ElementAx;
#include <ElementAx.h>

/** Axisymetric finite element class.
This class implements a 4 nodes quadrilateral axisymetric finite element.
see Domain, Node, IntegrationPoint, ElementAx
author Olivier PANTALE
version DynELA v. 0.9.3
*/
/** @dia:pos 86.8003,206.8 */
/** @dia:route ElementAx;v,96.9503,199,202.85,98.9503,206.8 */
class ElQua4nAx : public ElementAx
{

public:
  /** @dia:route 3,16;h,86.8003,207.5,83.7503,221.85,18.75,144.9,41.8 */
  static const ElementData elementData;

public:
  /* constructeurs */
  ElQua4nAx(Indice No_ = 1);
  ElQua4nAx(const ElQua4nAx &el);
  ~ElQua4nAx();

  Indice numberOfUnderIntegrationPoints()
  {
    return 1;
  }

  //  void computeInternalMatrices (Boolean reference=False);

  Real getLength(); /* calcule et renvoie la longueur caracteristique */
  Real getVolume();
  // Boolean getUnderIntegrPointCoords (Indice, Vec3D & coords, Real& weight) ;
  void getShapeFunctionAtPoint(Vector &N, const Vec3D &point) const;
  void getDerShapeFunctionAtPoint(Matrix &N, const Vec3D &point);
  void computeGlob2Loc();
  void glob2Loc(const Vec3D &point, Vec3D &local);
  //void getIntgtoNodes (Vector & N, const Vec3D & point) const;
};

#endif
