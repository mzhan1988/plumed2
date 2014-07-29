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
#include "RMSDBase.h"

namespace PLMD{

RMSDBase::RMSDBase( const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
SingleDomainRMSD(ro)
{
}

double RMSDBase::calculate( const std::vector<Vector>& pos, const bool& squared ){
  clearDerivatives(); 
  return calc( pos, squared );
}    

double RMSDBase::calc( const std::vector<Vector>& pos, const Pbc& pbc, const bool& squared ){
  plumed_dbg_assert( pos.size()==getNumberOfAtoms() );
  return calc( pos, squared );
}

// this is just a function that does not prevent not to have this feature but warns you if you try to use without a proper implementation
double RMSDBase::calc_DDistDRef(const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef){
  plumed_merror("You need to define calc_DDistDRef for this metrics if you want to use it ");
  return 0.;
}
double RMSDBase::calc_DDistDRef_Rot_DRotDPos(const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos){
  plumed_merror("You need to define calc_DDistDRef_Rot_DRotDPos for this metrics if you want to use it ");
  return 0.;
}
double RMSDBase::calc_DDistDRef_Rot_DRotDPos_DRotDRef(const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,Matrix<std::vector<Vector> > & DRotDRef ){
  plumed_merror("You need to define calc_DDistDRef_Rot_DRotDPos_DRotDRef for this metrics if you want to use it ");
  return 0.;
}



// this is here just not to make fail the others: just in case they are not
double RMSDBase::calculate_DDistDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef){
  return calc_DDistDRef( pos, squared , DDistDRef);
}


}
