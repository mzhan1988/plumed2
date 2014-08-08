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
#include "MetricRegister.h"
#include "RMSDBase.h"
#include "tools/Matrix.h"
#include "tools/Exception.h"
#include <iostream>
using namespace std;

namespace PLMD{

/// this is a class which is needed to share information across the various non-threadsafe routines
/// so that the public function of rmsd are threadsafe while the inner core can safely share information 
class RMSDCoreData
{
	private:
		bool alEqDis;
		bool distanceIsSquared;
		bool hasDistance;
		bool isInitialized;
		bool safe;
// use initialization reference assignment to speed up instead of copying and sucking out memory
// Reference coordinates
                const std::vector<Vector> &positions;
                const std::vector<Vector> &reference;
                // Weights for alignment
                const std::vector<double> &align;
                // Weights for deviation
                const std::vector<double> &displace;
// the needed stuff for distance and more (one could use eigenvecs components and eigenvals for some reason)
		double dist;
		std::vector<double> eigenvals;
		Matrix<double> eigenvecs;
		double rr00; //  sum of positions squared (needed for dist calc)
		double rr11; //  sum of reference squared (needed for dist calc)
// rotation derived from the eigenvector having the smallest eigenvalue
		Tensor rotation;
// derivative of the rotation only available when align!=displace
		Tensor drotation_drr01[3][3];
		Tensor ddist_drr01;
	        Tensor ddist_drotation;
// difference of components
		std::vector<Vector> d;
// geometric center of the running position and reference
 		Vector cpositions,creference;
	public:
		// the constructor (note: only references are passed, therefore is rather fast)
		// note:: this aligns the reference onto the positions
		RMSDCoreData(const std::vector<double> &a ,const std::vector<double> &d,const std::vector<Vector> &p, const std::vector<Vector> &r ):alEqDis(false),distanceIsSquared(false),hasDistance(false),isInitialized(false),safe(false),positions(p),reference(r),align(a),displace(d){};
		//  does the core calc : first thing to call after the constructor	
		void doCoreCalc(bool safe,bool alEqDis);
		// retrieve the distance if required after doCoreCalc 
		double getDistance(bool squared);
		// retrieve the derivative of the distance respect to the position
		std::vector<Vector> getDDistanceDPositions();
		// retrieve the derivative of the distance respect to the reference
		std::vector<Vector> getDDistanceDReference();
		// get aligned reference onto position
                std::vector<Vector> getAlignedReferenceToPositions();	
		// get aligned position onto reference
                std::vector<Vector> getAlignedPositionsToReference();	
		// get centered positions
                std::vector<Vector> getCenteredPositions();	
		// get centered reference
                std::vector<Vector> getCenteredReference();	
		// get rotation matrix (reference ->positions) 
		Tensor getRotationMatrixReferenceToPositions();
		// get rotation matrix (positions -> reference) 
		Tensor getRotationMatrixPositionsToReference();
		// get the derivative of the rotation matrix respect to positions
		// note that the this transformation overlap the  reference onto position
		// if inverseTransform=true then aligns the positions onto reference
		Matrix<std::vector<Vector> > getDRotationDPosition( bool inverseTransform=false );
		// get the derivative of the rotation matrix respect to reference 
		// note that the this transformation overlap the  reference onto position
		// if inverseTransform=true then aligns the positions onto reference
		Matrix<std::vector<Vector> >  getDRotationDReference(bool inverseTransform=false );
};

class OptimalRMSD : public RMSDBase {
private:
  bool fast;
public:
  OptimalRMSD(const ReferenceConfigurationOptions& ro);
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );
  double calc_DDistDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef); 
  double calc_DDistDRef_Rot_DRotDPos( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos); 
  double calc_DDistDRef_Rot_DRotDPos_DRotDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,Matrix<std::vector<Vector> > & DRotDRef ); 
  double calc_PCA( const std::vector<Vector>& pos, const bool& squared , Tensor & Rotation, std::vector<Vector>  & ddistdpos , Matrix<std::vector<Vector> > & DRotDPos,std::vector<Vector>  & alignedpositions, std::vector<Vector> & centeredpositions, std::vector<Vector> &centeredreference); 

  template <bool safe,bool alEqDis>
  double optimalAlignment(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,
                          bool squared=false);

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,std::vector<Vector> & DDistDRef,
                          bool squared=false);

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef_Rot_DRotDPos(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,std::vector<Vector> & DDistDRef, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,
                          bool squared=false);

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,std::vector<Vector> & DDistDRef, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos, Matrix<std::vector<Vector> > & DRotDRef,
                          bool squared=false);

  template <bool safe,bool alEqDis>
  double optimalAlignment_PCA(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              std::vector<Vector> & alignedpositions,
                              std::vector<Vector> & centeredpositions, 
                              std::vector<Vector> & centeredreference, 
			      Tensor & Rotation, 
                              std::vector<Vector> & ddistdpos, 
			      Matrix<std::vector<Vector> > & DRotDPos,
                              bool squared=false);

};

PLUMED_REGISTER_METRIC(OptimalRMSD,"OPTIMAL")

OptimalRMSD::OptimalRMSD(const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
RMSDBase(ro)
{
  fast=ro.usingFastOption();
}

void OptimalRMSD::read( const PDB& pdb ){
  readReference( pdb ); 
}

double OptimalRMSD::calc( const std::vector<Vector>& pos, const bool& squared ){
  if( fast ){
     if( getAlign()==getDisplace() ) return optimalAlignment<false,true>(getAlign(),getDisplace(),pos,squared); 
     return optimalAlignment<false,false>(getAlign(),getDisplace(),pos,squared);
  } else {
     if( getAlign()==getDisplace() ) return optimalAlignment<true,true>(getAlign(),getDisplace(),pos,squared);
     return optimalAlignment<true,false>(getAlign(),getDisplace(),pos,squared);
  }
}

double OptimalRMSD::calc_DDistDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef ){
  if( fast ){
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef<false,true>(getAlign(),getDisplace(),pos,DDistDRef,squared); 
     return optimalAlignment_DDistDRef<false,false>(getAlign(),getDisplace(),pos,DDistDRef,squared);
  } else {
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef<true,true>(getAlign(),getDisplace(),pos,DDistDRef,squared);
     return optimalAlignment_DDistDRef<true,false>(getAlign(),getDisplace(),pos,DDistDRef,squared);
  }
}

double OptimalRMSD::calc_DDistDRef_Rot_DRotDPos( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos ){
  if( fast ){
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef_Rot_DRotDPos<false,true>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,squared); 
     return optimalAlignment_DDistDRef_Rot_DRotDPos<false,false>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,squared);
  } else {
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef_Rot_DRotDPos<true,true>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,squared);
     return optimalAlignment_DDistDRef_Rot_DRotDPos<true,false>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,squared);
  }
}

double OptimalRMSD::calc_DDistDRef_Rot_DRotDPos_DRotDRef( const std::vector<Vector>& pos, const bool& squared , std::vector<Vector>& DDistDRef,Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos , Matrix<std::vector<Vector> > & DRotDRef){
  if( fast ){
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<false,true>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,DRotDRef,squared); 
     return optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<false,false>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,DRotDRef,squared);
  } else {
     if( getAlign()==getDisplace() ) return optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<true,true>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,DRotDRef,squared);
     return optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<true,false>(getAlign(),getDisplace(),pos,DDistDRef,Rotation,DRotDPos,DRotDRef,squared);
  }
}

double OptimalRMSD::calc_PCA( const std::vector<Vector>& pos, const bool& squared , Tensor & Rotation, std::vector<Vector>   & ddistdpos ,Matrix<std::vector<Vector> > & DRotDPos, std::vector<Vector>   & alignedpositions, std::vector<Vector>  & centeredpositions, std::vector<Vector>  & centeredreference){
//  if( fast ){
     //if( getAlign()==getDisplace() ) return optimalAlignment_PCA<false,true>(getAlign(),getDisplace(),pos,alignedpositions,centeredpositions,centeredreference,Rotation,DRotDPos,squared); 
//     return optimalAlignment_PCA<false,false>(getAlign(),getDisplace(),pos,alignedpositions, centeredpositions,centeredreference,Rotation,DRotDPos,squared);
//  } else {
     //if( getAlign()==getDisplace() ) return optimalAlignment_PCA<true,true>(getAlign(),getDisplace(),pos,alignedpositions,centeredpositions,centeredreference,Rotation,DRotDPos,squared);
     return optimalAlignment_PCA<true,false>(getAlign(),getDisplace(),pos,alignedpositions, centeredpositions,centeredreference,Rotation,ddistdpos,DRotDPos,squared);
///  }
}



// standard RMSD calculation, just calculate
#ifdef OLDRMSD
template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              bool squared){
  double dist(0);
  unsigned n=getNumberOfReferencePositions();
// This is the trace of positions*positions + reference*reference
  double rr00(0);
  double rr11(0);
// This is positions*reference
  Tensor rr01; Vector rpos;

  Vector cpositions;

// first expensive loop: compute centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    cpositions+=positions[iat]*w;
  }

// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat]; rpos=getReferencePosition(iat);
    rr00+=dotProduct(positions[iat]-cpositions,positions[iat]-cpositions)*w;
    rr11+=dotProduct(rpos,rpos)*w;
    rr01+=Tensor(positions[iat]-cpositions,rpos)*w;
  }

  Matrix<double> m=Matrix<double>(4,4);
  m[0][0]=2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
  m[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
  m[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
  m[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
  m[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
  m[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
  m[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];

  Tensor dm_drr01[4][4];
  if(!alEqDis){
    dm_drr01[0][0] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[1][1] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[2][2] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[3][3] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[0][1] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,+1.0, 0.0);
    dm_drr01[0][2] = 2.0*Tensor( 0.0, 0.0,+1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[0][3] = 2.0*Tensor( 0.0,-1.0, 0.0, +1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][2] = 2.0*Tensor( 0.0,-1.0, 0.0, -1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][3] = 2.0*Tensor( 0.0, 0.0,-1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[2][3] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,-1.0, 0.0);
    dm_drr01[1][0] = dm_drr01[0][1];
    dm_drr01[2][0] = dm_drr01[0][2];
    dm_drr01[2][1] = dm_drr01[1][2];
    dm_drr01[3][0] = dm_drr01[0][3];
    dm_drr01[3][1] = dm_drr01[1][3];
    dm_drr01[3][2] = dm_drr01[2][3];
  }

  std::vector<double> eigenvals;
  Matrix<double> eigenvecs;
  int diagerror=diagMat(m, eigenvals, eigenvecs );

  if (diagerror!=0){
    std::string sdiagerror;
    Tools::convert(diagerror,sdiagerror);
    std::string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
    plumed_merror(msg);
  }

  dist=eigenvals[0]+rr00+rr11;

  Matrix<double> ddist_dm(4,4);

  Vector4d q(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);

  Tensor dq_drr01[4];
  if(!alEqDis){
    double dq_dm[4][4][4];
    for(unsigned i=0;i<4;i++) for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++){
      double tmp=0.0;
// perturbation theory for matrix m
      for(unsigned l=1;l<4;l++) tmp+=eigenvecs[l][j]*eigenvecs[l][i]/(eigenvals[0]-eigenvals[l])*eigenvecs[0][k];
      dq_dm[i][j][k]=tmp;
    }
// propagation to _drr01
    for(unsigned i=0;i<4;i++){
      Tensor tmp;
      for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++) {
        tmp+=dq_dm[i][j][k]*dm_drr01[j][k];
      }
      dq_drr01[i]=tmp;
    }
  }

// This is the rotation matrix that brings reference to positions
// i.e. matmul(rotation,reference[iat])+shift is fitted to positions[iat]

  Tensor rotation;
  rotation[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotation[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotation[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  rotation[0][1]=2*(+q[0]*q[3]+q[1]*q[2]);
  rotation[0][2]=2*(-q[0]*q[2]+q[1]*q[3]);
  rotation[1][2]=2*(+q[0]*q[1]+q[2]*q[3]);
  rotation[1][0]=2*(-q[0]*q[3]+q[1]*q[2]);
  rotation[2][0]=2*(+q[0]*q[2]+q[1]*q[3]);
  rotation[2][1]=2*(-q[0]*q[1]+q[2]*q[3]);


  Tensor drotation_drr01[3][3];
  if(!alEqDis){
    drotation_drr01[0][0]=2*q[0]*dq_drr01[0]+2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[1][1]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]+2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[2][2]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]+2*q[3]*dq_drr01[3];
    drotation_drr01[0][1]=2*(+(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[0][2]=2*(-(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[1][2]=2*(+(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
    drotation_drr01[1][0]=2*(-(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[2][0]=2*(+(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[2][1]=2*(-(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
  }

  double prefactor=2.0;

  if(!squared && alEqDis) prefactor*=0.5/sqrt(dist);

// if "safe", recompute dist here to a better accuracy
  if(safe || !alEqDis) dist=0.0;

// If safe is set to "false", MSD is taken from the eigenvalue of the M matrix
// If safe is set to "true", MSD is recomputed from the rotational matrix
// For some reason, this last approach leads to less numerical noise but adds an overhead

  Tensor ddist_drotation;
  Vector ddist_dcpositions;

// third expensive loop: derivatives
  for(unsigned iat=0;iat<n;iat++){
    Vector d(positions[iat]-cpositions - matmul(rotation,getReferencePosition(iat)));
    if(alEqDis){
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
       addAtomicDerivatives( iat, prefactor*align[iat]*d );
       if(safe) dist+=align[iat]*modulo2(d);
    } else {
// the case for align != displace is different, sob:
      dist+=displace[iat]*modulo2(d);
// these are the derivatives assuming the roto-translation as frozen
      atom_ders[iat]=2*displace[iat]*d;
// here I accumulate derivatives wrt rotation matrix ..
      ddist_drotation+=-2*displace[iat]*extProduct(d,getReferencePosition(iat));
// .. and cpositions
      ddist_dcpositions+=-2*displace[iat]*d;
    }
  }

  if(!alEqDis){
    Tensor ddist_drr01;
    for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) ddist_drr01+=ddist_drotation[i][j]*drotation_drr01[i][j];
    for(unsigned iat=0;iat<n;iat++){
// this is propagating to positions.
// I am implicitly using the derivative of rr01 wrt positions here
      addAtomicDerivatives( iat, matmul(ddist_drr01,getReferencePosition(iat))*align[iat] );
      addAtomicDerivatives( iat, ddist_dcpositions*align[iat] );
    }
  }
  if(!squared){
    dist=sqrt(dist);
    if(!alEqDis){
      double xx=0.5/dist;
      for(unsigned iat=0;iat<atom_ders.size();iat++) atom_ders[iat]*=xx;
    }
  }

  return dist;
}
#else
// this is the standard version: no renewal of reference
template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              bool squared){
   //initialize the data into the structure
   RMSDCoreData cd(align,displace,positions,getReferencePositions()); 
   // Perform the diagonalization and all the needed stuff
   cd.doCoreCalc(safe,alEqDis); 
   // make the core calc distance
   double dist=cd.getDistance(squared); 
   // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...) 
   atom_ders=cd.getDDistanceDPositions(); 
   return dist;    
}

template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment_DDistDRef(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              std::vector<Vector> & ddistdref,
                              bool squared){
   //initialize the data into the structure
   RMSDCoreData cd(align,displace,positions,getReferencePositions()); 
   // Perform the diagonalization and all the needed stuff
   cd.doCoreCalc(safe,alEqDis); 
   // make the core calc distance
   double dist=cd.getDistance(squared); 
   // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...) 
   atom_ders=cd.getDDistanceDPositions(); 
   // get also the derivative respect to the reference 
   ddistdref=cd.getDDistanceDReference(); 
   return dist;    
}


template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment_DDistDRef_Rot_DRotDPos(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              std::vector<Vector> & ddistdref, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,
                              bool squared){
   //initialize the data into the structure
   RMSDCoreData cd(align,displace,positions,getReferencePositions()); 
   // Perform the diagonalization and all the needed stuff
   cd.doCoreCalc(safe,alEqDis); 
   // make the core calc distance
   double dist=cd.getDistance(squared); 
   // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...) 
   atom_ders=cd.getDDistanceDPositions(); 
   // get also the derivative respect to the reference 
   ddistdref=cd.getDDistanceDReference(); 
   // get the rotation matrix
   Rotation=cd.getRotationMatrixReferenceToPositions(); 
   // get its derivative
   DRotDPos=cd.getDRotationDPosition();  
   return dist;    
}

template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              std::vector<Vector> & ddistdref, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,Matrix<std::vector<Vector> > & DRotDRef,
                              bool squared){
   //initialize the data into the structure
   RMSDCoreData cd(align,displace,positions,getReferencePositions()); 
   // Perform the diagonalization and all the needed stuff
   cd.doCoreCalc(safe,alEqDis); 
   // make the core calc distance
   double dist=cd.getDistance(squared); 
   // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...) 
   atom_ders=cd.getDDistanceDPositions(); 
   // get also the derivative respect to the reference 
   ddistdref=cd.getDDistanceDReference(); 
   // get the rotation matrix
   Rotation=cd.getRotationMatrixReferenceToPositions(); 
   // get its derivative
   DRotDPos=cd.getDRotationDPosition();  
   DRotDRef=cd.getDRotationDReference();  
   return dist;    
}


/*
This below is a convenience call to calculate PCA components which are needed
Remember that in this case the Reference is aligned onto the positions
therefore the pca should read ( x_t-com_t -R (x_r-com_r) ) * Rv
*/
template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment_PCA(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              std::vector<Vector> & alignedpositions,
                              std::vector<Vector> & centeredpositions, 
                              std::vector<Vector> & centeredreference, 
			      Tensor & Rotation, 
                              std::vector<Vector> & ddistdpos,
			      Matrix<std::vector<Vector> > & DRotDPos,
                              bool squared){
   //initialize the data into the structure
   RMSDCoreData cd(align,displace,positions,getReferencePositions()); 
   // Perform the diagonalization and all the needed stuff
   cd.doCoreCalc(safe,alEqDis); 
   // make the core calc distance
   double dist=cd.getDistance(squared); 
   // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
   ddistdpos=cd.getDDistanceDPositions(); 	
   // get the rotation matrix
   Rotation=cd.getRotationMatrixPositionsToReference(); 
   // get its derivative
   DRotDPos=cd.getDRotationDPosition(true); // this gives back the inverse  
   // get aligned positions 
   alignedpositions=cd.getAlignedPositionsToReference(); 
   // get centered positions
   centeredpositions=cd.getCenteredPositions();
   // get centered reference
   centeredreference=cd.getCenteredReference();

   return dist;    
}



/// This calculates the elements needed by the quaternion to calculate everything that is needed
/// additional calls retrieve different components
void RMSDCoreData::doCoreCalc(bool safe,bool alEqDis){

  const unsigned n=static_cast<unsigned int>(reference.size());

  cpositions.zero();
  creference.zero();

// first expensive loop: compute centers
// TODO: additional flags could avoid it in case the reference is already centered 
// or create a specific object for reference and position that calculates the com
// when the object is created and recalculates it whenever is touched by some public method
// (this last one is not so smart since then one must ensure that in the ref and pos the weights are identical 
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    cpositions+=positions[iat]*w;
    creference+=reference[iat]*w;
  }

// This is the trace of positions*positions + reference*reference
  rr00=0.;
  rr11=0.;
// This is positions*reference
  Tensor rr01;
// second expensive loop: compute second moments wrt centers
// TODO is  w=align[iat] actually right???
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    rr00+=dotProduct(positions[iat]-cpositions,positions[iat]-cpositions)*w;
    rr11+=dotProduct(reference[iat]-creference,reference[iat]-creference)*w;
    rr01+=Tensor(positions[iat]-cpositions,reference[iat]-creference)*w;
  }

// the quaternion matrix: this is internal
  Matrix<double> m=Matrix<double>(4,4);

  m[0][0]=2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
  m[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
  m[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
  m[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
  m[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
  m[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
  m[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];

  
  Tensor dm_drr01[4][4];
  if(!alEqDis){
    dm_drr01[0][0] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,-1.0); 
    dm_drr01[1][1] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[2][2] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[3][3] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[0][1] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,+1.0, 0.0);
    dm_drr01[0][2] = 2.0*Tensor( 0.0, 0.0,+1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[0][3] = 2.0*Tensor( 0.0,-1.0, 0.0, +1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][2] = 2.0*Tensor( 0.0,-1.0, 0.0, -1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][3] = 2.0*Tensor( 0.0, 0.0,-1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[2][3] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,-1.0, 0.0);
    dm_drr01[1][0] = dm_drr01[0][1];
    dm_drr01[2][0] = dm_drr01[0][2];
    dm_drr01[2][1] = dm_drr01[1][2];
    dm_drr01[3][0] = dm_drr01[0][3];
    dm_drr01[3][1] = dm_drr01[1][3];
    dm_drr01[3][2] = dm_drr01[2][3];
  }


  int diagerror=diagMat(m, eigenvals, eigenvecs );

  if (diagerror!=0){
    std::string sdiagerror;
    Tools::convert(diagerror,sdiagerror);
    std::string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
    plumed_merror(msg);
  }

  Vector4d q(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);

  Tensor dq_drr01[4];
  if(!alEqDis){
    double dq_dm[4][4][4];
    for(unsigned i=0;i<4;i++) for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++){
      double tmp=0.0;
// perturbation theory for matrix m
      for(unsigned l=1;l<4;l++) tmp+=eigenvecs[l][j]*eigenvecs[l][i]/(eigenvals[0]-eigenvals[l])*eigenvecs[0][k];
      dq_dm[i][j][k]=tmp;
    }
// propagation to _drr01
    for(unsigned i=0;i<4;i++){
      Tensor tmp;
      for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++) {
        tmp+=dq_dm[i][j][k]*dm_drr01[j][k];
      }
      dq_drr01[i]=tmp;
    }
  }

// This is the rotation matrix that brings reference to positions
// i.e. matmul(rotation,reference[iat])+shift is fitted to positions[iat]

  rotation[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotation[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotation[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  rotation[0][1]=2*(+q[0]*q[3]+q[1]*q[2]);
  rotation[0][2]=2*(-q[0]*q[2]+q[1]*q[3]);
  rotation[1][2]=2*(+q[0]*q[1]+q[2]*q[3]);
  rotation[1][0]=2*(-q[0]*q[3]+q[1]*q[2]);
  rotation[2][0]=2*(+q[0]*q[2]+q[1]*q[3]);
  rotation[2][1]=2*(-q[0]*q[1]+q[2]*q[3]);
  
  if(!alEqDis){
    drotation_drr01[0][0]=2*q[0]*dq_drr01[0]+2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[1][1]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]+2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[2][2]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]+2*q[3]*dq_drr01[3];
    drotation_drr01[0][1]=2*(+(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[0][2]=2*(-(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[1][2]=2*(+(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
    drotation_drr01[1][0]=2*(-(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[2][0]=2*(+(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[2][1]=2*(-(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
  }

  d.resize(n);

  // calculate rotation matrix derivatives and components distances needed for components only when align!=displacement
  if(!alEqDis)ddist_drotation.zero();
  for(unsigned iat=0;iat<n;iat++){
    // components differences: this is useful externally
    d[iat]=positions[iat]-cpositions - matmul(rotation,reference[iat]-creference);	
    // ddist_drotation if needed
    if(!alEqDis) ddist_drotation+=-2*displace[iat]*extProduct(d[iat],reference[iat]-creference);
  }

  if(!alEqDis){
          ddist_drr01.zero();
          for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) ddist_drr01+=ddist_drotation[i][j]*drotation_drr01[i][j];
  }
  // transfer this bools to the cd so that this settings will be reflected in the other calls
  this->alEqDis=alEqDis; 
  this->safe=safe; 
  isInitialized=true;

}
/// just retrieve the distance already calculated
double RMSDCoreData::getDistance(bool squared){

  if(!isInitialized)plumed_merror("OptimalRMSD.cpp cannot calculate the distance without being initialized first by doCoreCalc ");
  dist=eigenvals[0]+rr00+rr11;
  double prefactor=2.0;
  
  if(!squared && alEqDis) prefactor*=0.5/sqrt(dist);
  if(safe || !alEqDis) dist=0.0;
  const unsigned n=static_cast<unsigned int>(reference.size());
  for(unsigned iat=0;iat<n;iat++){
  	if(alEqDis){
  	    if(safe) dist+=align[iat]*modulo2(d[iat]);
  	} else {
  	    dist+=displace[iat]*modulo2(d[iat]);
  	}
  }
  if(!squared){
  	dist=sqrt(dist);
	distanceIsSquared=false;
  }else{
	distanceIsSquared=true;
  }
  hasDistance=true;
  return dist; 
}

std::vector<Vector> RMSDCoreData::getDDistanceDPositions(){
  std::vector<Vector>  derivatives;
  const unsigned n=static_cast<unsigned int>(reference.size());
  Vector ddist_dcpositions;
  derivatives.resize(n);
  double prefactor=2.0;
  if(!hasDistance)plumed_merror("getDPositionsDerivatives needs to calculate the distance via getDistance first !");
  if(!isInitialized)plumed_merror("getDPositionsDerivatives needs to initialize the coreData first!");
  vector<Vector> ddist_tmp(n);
  Vector csum;
  Vector tmp1,tmp2;
  for(unsigned iat=0;iat<n;iat++){
    if(alEqDis){
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      derivatives[iat]= prefactor*align[iat]*d[iat];
    } else {
// these are the derivatives assuming the roto-translation as frozen
      tmp1=2*displace[iat]*d[iat];
      derivatives[iat]=tmp1;
// derivative of cpositions
      ddist_dcpositions+=-tmp1;
      // these needed for com corrections
      tmp2=matmul(ddist_drr01,reference[iat]-creference)*align[iat];	
      derivatives[iat]+=tmp2;	
      csum+=tmp2;
    }
  }

  if(!alEqDis){
    for(unsigned iat=0;iat<n;iat++)derivatives[iat]+=(ddist_dcpositions-csum)*align[iat]; 
  }
  if(!distanceIsSquared){
    if(!alEqDis){
      double xx=0.5/dist;
      for(unsigned iat=0;iat<n;iat++) derivatives[iat]*=xx;
    }
  }
  return derivatives;
}

std::vector<Vector>  RMSDCoreData::getDDistanceDReference(){
  std::vector<Vector>  derivatives;
  const unsigned n=static_cast<unsigned int>(reference.size());
  Vector ddist_dcreference;
  derivatives.resize(n);
  double prefactor=2.0;
  vector<Vector> ddist_tmp(n);
  Vector csum,tmp1,tmp2;

  if(!hasDistance)plumed_merror("getDDistanceDReference needs to calculate the distance via getDistance first !");
  if(!isInitialized)plumed_merror("getDDistanceDReference to initialize the coreData first!");
  // get the transpose rotation
  Tensor t_rotation=rotation.transpose();
  Tensor t_ddist_drr01=ddist_drr01.transpose();	
  
// third expensive loop: derivatives
  for(unsigned iat=0;iat<n;iat++){
    if(alEqDis){
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
	//TODO: check this derivative down here
      derivatives[iat]= -prefactor*align[iat]*matmul(t_rotation,d[iat]);
    } else {
// these are the derivatives assuming the roto-translation as frozen
      tmp1=2*displace[iat]*matmul(t_rotation,d[iat]);
      derivatives[iat]= -tmp1;
// derivative of cpositions
      ddist_dcreference+=tmp1;
      // these below are needed for com correction
      tmp2=matmul(t_ddist_drr01,positions[iat]-cpositions)*align[iat];
      derivatives[iat]+=tmp2;	
      csum+=tmp2; 
    }
  }

  if(!alEqDis){
    for(unsigned iat=0;iat<n;iat++)derivatives[iat]+=(ddist_dcreference-csum)*align[iat]; 
  }
  if(!distanceIsSquared){
    if(!alEqDis){
      double xx=0.5/dist;
      for(unsigned iat=0;iat<n;iat++) derivatives[iat]*=xx;
    }
  }
  return derivatives;
}

/*
This below is the derivative of the rotation matrix that aligns the reference onto the positions
respect to positions
note that the this transformation overlap the  reference onto position
if inverseTransform=true then aligns the positions onto reference
*/
Matrix<std::vector<Vector> >  RMSDCoreData::getDRotationDPosition( bool inverseTransform ){
  const unsigned n=static_cast<unsigned int>(reference.size());
  if(!isInitialized)plumed_merror("getDRotationDPosition to initialize the coreData first!");
  Matrix<std::vector<Vector> > DRotDPos=Matrix<std::vector<Vector> >(3,3);  
  // remember drotation_drr01 is Tensor drotation_drr01[3][3]
  //           (3x3 rot) (3x3 components of rr01)    
  std::vector<Vector> v(n);
  Vector csum;
  // these below could probably be calculated in the main routine
  for(unsigned iat=0;iat<n;iat++) csum+=(reference[iat]-creference)*align[iat];
  for(unsigned iat=0;iat<n;iat++) v[iat]=(reference[iat]-creference-csum)*align[iat];
  for(unsigned a=0;a<3;a++){
  	for(unsigned b=0;b<3;b++){
		if(inverseTransform){
			DRotDPos[b][a].resize(n);
			for(unsigned iat=0;iat<n;iat++){
  			      DRotDPos[b][a][iat]=matmul(drotation_drr01[a][b],v[iat]);
  			}
		}else{
			DRotDPos[a][b].resize(n);
  			for(unsigned iat=0;iat<n;iat++){
  			      DRotDPos[a][b][iat]=matmul(drotation_drr01[a][b],v[iat]);
  			}
  		}
  	}
  }
  return DRotDPos;
}

/*
This below is the derivative of the rotation matrix that aligns the reference onto the positions
respect to reference
note that the this transformation overlap the  reference onto position
if inverseTransform=true then aligns the positions onto reference
*/
Matrix<std::vector<Vector> >  RMSDCoreData::getDRotationDReference( bool inverseTransform ){
  const unsigned n=static_cast<unsigned int>(reference.size());
  if(!isInitialized)plumed_merror("getDRotationDPosition to initialize the coreData first!");
  Matrix<std::vector<Vector> > DRotDRef=Matrix<std::vector<Vector> >(3,3);  
  // remember drotation_drr01 is Tensor drotation_drr01[3][3]
  //           (3x3 rot) (3x3 components of rr01)    
  std::vector<Vector> v(n);
  Vector csum;
  // these below could probably be calculated in the main routine
  for(unsigned iat=0;iat<n;iat++) csum+=(positions[iat]-cpositions)*align[iat];
  for(unsigned iat=0;iat<n;iat++) v[iat]=(positions[iat]-cpositions-csum)*align[iat];
 
  for(unsigned a=0;a<3;a++){
  	for(unsigned b=0;b<3;b++){
		Tensor t_drotation_drr01=drotation_drr01[a][b].transpose();
		if(inverseTransform){
			DRotDRef[b][a].resize(n);
			for(unsigned iat=0;iat<n;iat++){
  			      DRotDRef[b][a][iat]=matmul(t_drotation_drr01,v[iat]);
  			}
		}else{
			DRotDRef[a][b].resize(n);
  			for(unsigned iat=0;iat<n;iat++){
  			      DRotDRef[a][b][iat]=matmul(t_drotation_drr01,v[iat]);
  			}
  		}
  	}
  }
  return DRotDRef;
}


std::vector<Vector> RMSDCoreData::getAlignedReferenceToPositions(){
	  std::vector<Vector> alignedref;
	  const unsigned n=static_cast<unsigned int>(reference.size());
	  alignedref.resize(n);
	  if(!isInitialized)plumed_merror("getAlignedReferenceToPostions needs to initialize the coreData first!");
          // avoid to calculate matrix element but use the sum of what you have		  
	  for(unsigned iat=0;iat<n;iat++)alignedref[iat]=-d[iat]+positions[iat]-cpositions;
	  return alignedref; 
}
std::vector<Vector> RMSDCoreData::getAlignedPositionsToReference(){
	  std::vector<Vector> alignedpos;
	  const unsigned n=static_cast<unsigned int>(positions.size());
	  alignedpos.resize(n);
	  if(!isInitialized)plumed_merror("getAlignedPostionsToReference needs to initialize the coreData first!");
          // avoid to calculate matrix element but use the sum of what you have		  
	  for(unsigned iat=0;iat<n;iat++)alignedpos[iat]=matmul(rotation.transpose(),positions[iat]-cpositions);
	  return alignedpos; 
}


std::vector<Vector> RMSDCoreData::getCenteredPositions(){
	  std::vector<Vector> centeredpos;
	  const unsigned n=static_cast<unsigned int>(reference.size());
	  centeredpos.resize(n);
	  if(!isInitialized)plumed_merror("getCenteredPositions needs to initialize the coreData first!");
          // avoid to calculate matrix element but use the sum of what you have		  
	  for(unsigned iat=0;iat<n;iat++)centeredpos[iat]=positions[iat]-cpositions;
	  return centeredpos; 
}

std::vector<Vector> RMSDCoreData::getCenteredReference(){
	  std::vector<Vector> centeredref;
	  const unsigned n=static_cast<unsigned int>(reference.size());
	  centeredref.resize(n);
	  if(!isInitialized)plumed_merror("getCenteredReference needs to initialize the coreData first!");
          // avoid to calculate matrix element but use the sum of what you have		  
	  for(unsigned iat=0;iat<n;iat++)centeredref[iat]=reference[iat]-creference;
	  return centeredref; 
}




Tensor RMSDCoreData::getRotationMatrixReferenceToPositions(){
	  if(!isInitialized)plumed_merror("getRotationMatrixReferenceToPositions needs to initialize the coreData first!");
	  return rotation;
}

Tensor RMSDCoreData::getRotationMatrixPositionsToReference(){
	  if(!isInitialized)plumed_merror("getRotationMatrixReferenceToPositions needs to initialize the coreData first!");
	  return rotation.transpose();
}




#endif


}
