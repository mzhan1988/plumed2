/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_reference_RMSDBase_h
#define __PLUMED_reference_RMSDBase_h

#include <vector>
#include <string>
#include "SingleDomainRMSD.h"
#include "tools/Matrix.h"

namespace PLMD {

class Pbc;

class RMSDBase : public SingleDomainRMSD {
public:
  RMSDBase( const ReferenceConfigurationOptions& ro );
  double calc( const std::vector<Vector>& pos, const Pbc& pbc, const bool& squared );
  double calculate( const std::vector<Vector>& pos, const bool& squared );
  double calculate_DDistDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef);
  virtual double calc_DDistDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef);
  virtual double calc_DDistDRef_Rot_DRotDPos( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos);
  virtual double calc_DDistDRef_Rot_DRotDPos_DRotDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos, Matrix<std::vector<Vector> > & DRotDRef);
  virtual double calc_PCA(const std::vector<Vector>& pos, const bool& squared , Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos, std::vector<Vector>   & alignedpositions, std::vector<Vector>  & centeredpositions ) ;
  virtual double calc( const std::vector<Vector>& pos, const bool& squared )=0;
};

}
#endif
