#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "plumed/tools/Exception.h"
#include "plumed/tools/RMSD.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/PDB.h"
#include "plumed/tools/Log.h"
#include "plumed/reference/RMSDBase.h"
#include "plumed/reference/MetricRegister.h"
#include "plumed/tools/Vector.h"
#include "plumed/tools/Matrix.h"
#include <vector>
#include <fstream>
#include "plumed/tools/AtomNumber.h"


using namespace std ;
using namespace PLMD;

int main(int argc, char* argv[]) {

  // findiff step: eps
  double eps=1.e-6;

  // default task, calculate the RMSD of one frame respect to a set of others 
  vector<int> task(1);task[0]=0;
  // first parse the task
  for(int i = 1; i < argc; i++){ 
      task.push_back(atoi(argv[i]));
  }  
  if(argc==1){
	cout<<"ARGUMENTS: \n";
        cout<<" -1 : squared=true (default false)\n";  
	cout<<" -2 : normalize_weights=false (default true)\n";
	cout<<" -3 : reset_com=false (default true) \n";
	cout<<"  0 : normal rmsd/msd calculation  and derivative dumps (default: always done)\n";
	cout<<"  1 : findiff test for  d msd / d position  (inhomogenehous weights)\n";
	cout<<"  2 : findiff test for  d msd / d reference (inhomogenehous weights)\n";
	cout<<"  3 : findiff test for  d msd / d position  (homogenehous weights)\n";
	cout<<"  4 : findiff test for  d msd / d reference (homogenehous weights)\n";
	cout<<"  5 : findiff test for  d Rot / d position  (inhomogenehous weights)  (reference->position)\n";
	cout<<"  6 : findiff test for  d Rot / d reference  (inhomogenehous weights) (reference->position)\n";
	cout<<"  7 : consistency check for MSD proportionality (works with squared=true through option -1 )\n";
	cout<<"  8 : do some timings for all the above routines and for a growing number of atoms\n";
	cout<<"  9 : test the rotation order: print position.pdb reference.pdb aligned.pdb and check that it makes sense(should be reference aligned onto positions)\n";
	cout<<" 10 : findiff test for  d Rot / d position  ( position -> reference ) \n";

	//return 0 ;
  }


  bool reset_com=true;
  bool normalize_weights=true;
  bool squared=false;

  if(std::find(task.begin(), task.end(), -1)!=task.end()){cout<<"squared=true (default false)"<<endl;squared=true;}
  if(std::find(task.begin(), task.end(), -2)!=task.end()){cout<<"normalize_weights=false (default true)"<<endl;normalize_weights=false;}
  if(std::find(task.begin(), task.end(), -3)!=task.end()){cout<<"reset_com=false (default true)"<<endl; reset_com=false;}

  PDB pdbref; 

  // now create the object
  PLMD::RMSDBase* rmsd;
  string type; type.assign("OPTIMAL");

  string reference; reference.assign("1GB1_mdl1_rototranslated.pdb");
  PDB pdb;

  // mimick gromacs: do it in nm
  if( !pdb.read(reference,false,0.1) ) 
        cout<<"missing input file 1GB1_mdl1_rototranslated.pdb "<<"\n";
 
  cout<<"NOW CREATING THE RMSD OBJECT...";  

  // for my reference: 
  // the register contains all the constructors: in this case 
  // create calls the constructor for "OPTIMAL" which is OptimalRMSD 
  // the templates defines which is the return type. In this case RMSDBase which is responsible to 
  // set the reference configuration (in this case calls RMSDBase::public SingleDomainRMSD 
  //                                                                       SingleDomainRMSD : public ReferenceAtoms 
  //                                                                       ReferenceAtoms : public ReferenceConfiguration
  // which has ReferenceConfiguration::set( const PDB& )
  // that calls a virtual  read( const PDB& )=0;  -> the implementation is in void OptimalRMSD::read( const PDB& pdb ) 
  // which calls     SingleDomainRMSD::readReference( const PDB& pdb );
  // this reads the atoms and reset to reference
  //
  rmsd = metricRegister().create<RMSDBase>(type,pdb);
  std::vector<AtomNumber> atoms;
  rmsd->getAtomRequests( atoms );
  rmsd->setNumberOfAtoms( atoms.size() );
  cout<<"DONE!"<<endl;  

  // now calculate something: take another conformation
  PDB pdbrun; 

  // mimic gromacs: do it in nm
  if( !pdbrun.read("1GB1_mdl2.pdb",false,0.1) )
        cout<<"missing input file 1GB1_mdl2.pdb\n" ;

  std::vector<Vector> run ;  run=pdbrun.getPositions() ;
  std::vector<Vector> ref ;  ref=pdb.getPositions() ;
  std::vector<double> align; align=pdb.getOccupancy();
  std::vector<double> displace;  displace=pdb.getBeta();

  // reset some stuff

  rmsd->setResetCom(reset_com);
  rmsd->setNormalizeWeights(normalize_weights);
  rmsd->setReferenceAtoms( ref, align, displace );



  
  // Task 0: calculate the alignment and dump some data
  if(std::find(task.begin(), task.end(), 0)!=task.end()){
	cout<<"Task 0: calculates the alignment and retrieve some data"<<endl;
  	double r=rmsd->calculate( run, squared ); 
	cout<<"RMSD IS "<<r<<"\n";	
	// now dump some more information
        ofstream myfile;
        myfile.open ("output_rmsd");
	myfile<<"RMSD ";
        myfile.setf( std::ios::fixed, std:: ios::floatfield );
	myfile.width(12);myfile.precision(6); myfile<<std::right<<r<<"\n";
	// the derivatives
	for(unsigned int i=0;i<run.size();++i){
	 myfile<<"DDIST_DPOS "<<std::setw(12)<<std::right<<rmsd->getAtomDerivative(i)[0]<<" "<<std::setw(12)<<std::right<<rmsd->getAtomDerivative(i)[1]<<" "<<std::setw(12)<<std::right<<rmsd->getAtomDerivative(i)[2]<<"\n";
	}
        std::vector<Vector> DDistDRef;	
	r=rmsd->calculate_DDistDRef( run, squared ,DDistDRef); 
	for(unsigned int i=0;i<run.size();++i){
	 myfile<<"DDIST_DREF "<<std::setw(12)<<std::right<<DDistDRef[i][0]<<" "<<std::setw(12)<<std::right<<DDistDRef[i][1]<<" "<<std::setw(12)<<std::right<<DDistDRef[i][2]<<"\n";
	}
	
	myfile.close();

  }
  // Task 1: calculate findiff of running frame
  if(std::find(task.begin(), task.end(), 1)!=task.end()){
	cout<<"Task 1: calculates the finite difference for the running frame"<<endl;
  	double r_old=rmsd->calculate( run, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			run[i][comp]+=eps;
  			double r=rmsd->calculate( run, squared ); 
			cout<<"DDIST_DPOS COMPONENT "<<comp<<" "<<(r-r_old)/(run[i][comp]-run_save[i][comp])<<" "<<rmsd->getAtomDerivative(i)[comp]<<"\n";
			// restore the old position
			run=run_save;
		}
	}
  }
  // Task 2: calculate findiff of reference frame
  if(std::find(task.begin(), task.end(), 2)!=task.end()){
	cout<<"Task 2: calculates the finite difference for the reference frame"<<endl;
  	double r_old=rmsd->calculate( run, squared ); 
	std::vector<Vector> ref_save=ref;	
	std::vector<Vector> DDistDRef;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			ref[i][comp]+=eps;
			// this function below also reset the com of the reference (but not of the running frame)
			rmsd->setReferenceAtoms( ref, align, displace );
  			double r=rmsd->calculate_DDistDRef( run, squared ,DDistDRef); 
			cout<<"DDIST_DREF COMPONENT "<<comp<<" "<<(r-r_old)/(ref[i][comp]-ref_save[i][comp])<<" "<<DDistDRef[i][comp]<<"\n";
			// restore the old position
			ref=ref_save;
		}
	}
  }
  // Task 3 calculate findiff of running frame for fast version (align=displace)
  std::vector<double> newalign(run.size(),1.); 
  std::vector<double> newdisplace(run.size(),1.); 
  rmsd->setReferenceAtoms( ref, newalign, newdisplace );
  if(std::find(task.begin(), task.end(), 3)!=task.end()){
	cout<<"Task 3: calculates the finite difference for the running frame (fast version)"<<endl;
  	double r_old=rmsd->calculate( run, squared ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			run[i][comp]+=eps;
  			double r=rmsd->calculate( run, squared ); 
			cout<<"DDIST_DPOS_FAST COMPONENT "<<comp<<" "<<(r-r_old)/(run[i][comp]-run_save[i][comp])<<" "<<rmsd->getAtomDerivative(i)[comp]<<"\n";
			// restore the old position
			run=run_save;
		}
	}
  }
  // Task 4: calculate findiff of reference frame for fast version (align=displace)
  if(std::find(task.begin(), task.end(), 4)!=task.end()){
	cout<<"Task 4: calculates the finite difference for the reference frame"<<endl;
	rmsd->setReferenceAtoms( ref, align, displace );
  	double r_old=rmsd->calculate( run, squared ); 
  	//r_old=rmsd->calculate( run, squared ); 
	//cout<<"ROLD "<<r_old<<endl;
	std::vector<Vector> DDistDRef;	
	std::vector<Vector> ref_save=ref;	
	for(unsigned int comp=0;comp<3;comp++){	
		for(unsigned int i=0;i<run.size();++i){
			//change the position		
			ref[i][comp]+=eps;
			// this function below also reset the com of the reference (but not of the running frame)
			rmsd->setReferenceAtoms( ref, align, displace );
  			double r=rmsd->calculate_DDistDRef( run, squared ,DDistDRef); 
			cout<<"DDIST_DREF_FAST COMPONENT "<<comp<<" "<<(r-r_old)/(ref[i][comp]-ref_save[i][comp])<<" "<<DDistDRef[i][comp]<<" "<<r<<" "<<r_old<<"\n";
			// restore the old position
			ref=ref_save;
		}
	}
  }
  // Task 5: calculate findiff of derivative of the rotation matrix respect to running fram
  rmsd->setReferenceAtoms( ref, align, displace );
  if(std::find(task.begin(), task.end(), 5)!=task.end()){
	cout<<"Task 5: calculates the finite difference for derivative of the rotation matrix respect to the the running frame"<<endl;
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
	rmsd->calc_DDistDRef_Rot_DRotDPos( run, squared, DDistDRef, OldRotation , DRotDPos  ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					run[i][comp]+=eps;
					rmsd->calc_DDistDRef_Rot_DRotDPos( run, squared, DDistDRef, Rotation , DRotDPos  ); 
					cout<<"DROT_DPOS COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(run[i][comp]-run_save[i][comp])<<" "<<DRotDPos[a][b][i][comp]<<"\n";
					// restore the old position
					run=run_save;
				}
			}
		}
	}
  }

  // Task 6: calculate findiff of derivative of the rotation matrix respect to reference frame 
  if(std::find(task.begin(), task.end(), 6)!=task.end()){
	cout<<"Task 6: calculates the finite difference for derivative of the rotation matrix respect to the the reference frame"<<endl;
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3),DRotDRef(3,3);
        std::vector<Vector> DDistDRef;
	rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, squared, DDistDRef, OldRotation , DRotDPos , DRotDRef ); 
	std::vector<Vector> ref_save=ref;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					ref[i][comp]+=eps;
					// this function below also reset the com of the reference (but not of the running frame)
					rmsd->setReferenceAtoms( ref, align, displace );
					rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( run, squared, DDistDRef, Rotation , DRotDPos, DRotDRef  ); 
					cout<<"DROT_DREF COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(ref[i][comp]-ref_save[i][comp])<<" "<<DRotDRef[a][b][i][comp]<<"\n";
					// restore the old position
					ref=ref_save;
				}
			}
		}
	}
  }
 
  // Task 7:  check weight consistency 

  if(std::find(task.begin(), task.end(), 7)!=task.end()){
	cout<<"Task 7: calculates the weight (displacement) consistency: all these should same result when weights are normalized in input by setReferenceAtoms otherwise they should be proportional when squared=true "<<endl;
  	double r=rmsd->calculate( run, squared ); 
	cout<<"STANDARD WEIGHT "<<r<<"\n"; 

        std::vector<double> newalign=align;//for(std::vector<double>::iterator p=newalign.begin();p!=newalign.end();++p ){(*p)*=2.;}  
        std::vector<double> newdisplace=displace;for(std::vector<double>::iterator p=newdisplace.begin();p!=newdisplace.end();++p ){(*p)*=2.;} 
	rmsd->setReferenceAtoms( ref, newalign, newdisplace );
  	r=rmsd->calculate( run, squared ); 
	cout<<"DOUBLE WEIGHT "<<r<<"\n"; 

        newalign=align;//for(std::vector<double>::iterator p=newalign.begin();p!=newalign.end();++p ){(*p)*=4.;}  
        newdisplace=displace;for(std::vector<double>::iterator p=newdisplace.begin();p!=newdisplace.end();++p ){(*p)*=4.;} 
	rmsd->setReferenceAtoms( ref, newalign, newdisplace );
  	r=rmsd->calculate( run, squared ); 
	cout<<"FOUR WEIGHT "<<r<<"\n"; 
  }

  // Task 8: do some timings to get a flavor
  if(std::find(task.begin(), task.end(), 8)!=task.end()){
	cout<<"Task 8: makes some timings for increasing atoms and different routines "<<endl;
	vector<Vector> r_run,r_ref;
	vector<double> r_al,r_disp;
	for (unsigned int i=0;i<10;i++){r_run.push_back(run[i]);r_ref.push_back(ref[i]);r_al.push_back(align[i]);r_disp.push_back(displace[i]);}

	for(unsigned int i=0;i<10;i++){
		cout<<"NUMBER OF ATOMS : "<<r_run.size()<<endl;
		unsigned ntest; ntest=100;
		// test the fast routine
	        rmsd->setReferenceAtoms( r_ref, r_al, r_disp );
		Stopwatch sw;
		sw.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calculate( r_run, squared );
		sw.stop();	
		cout<<"SIMPLE ROUTINE \n"<<sw<<endl;

	        std::vector<Vector> DDistDRef;
		Stopwatch sw2;
		sw2.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calculate_DDistDRef( r_run, squared ,DDistDRef); 
		sw2.stop();	
		cout<<"WITH REFERENCE FRAME: \n"<<sw2<<endl;

		Tensor Rotation;
        	Matrix<std::vector<Vector> > DRotDPos(3,3);	
		Stopwatch sw3;
		sw3.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calc_DDistDRef_Rot_DRotDPos( r_run, squared ,DDistDRef, Rotation , DRotDPos); 
		sw3.stop();	
		cout<<"WITH ROTATION MATRIX DERIVATIVE: \n"<<sw3<<endl;

                Matrix<std::vector<Vector> > DRotDRef(3,3);
		Stopwatch sw4;
		sw4.start();	
	        for(unsigned int j=0;j<ntest;j++)rmsd->calc_DDistDRef_Rot_DRotDPos_DRotDRef( r_run, squared ,DDistDRef, Rotation , DRotDPos, DRotDRef); 
		sw4.stop();	
		cout<<"WITH ROTATION MATRIX DERIVATIVE OF REEFERENCE: \n"<<sw4<<endl;
		// duplicate the atoms
		unsigned s=r_run.size();
		for (unsigned int i=0;i<s;i++){r_run.push_back(r_run[i]);r_ref.push_back(r_ref[i]);r_al.push_back(r_al[i]);r_disp.push_back(r_disp[i]);}
	
	}

  } 

  // Task 9: check the rotation
  if(std::find(task.begin(), task.end(), 9)!=task.end()){
	// dump the reference
	ofstream myfile;
	myfile.open ("reference.pdb");		
	std::vector<AtomNumber> at=pdb.getAtomNumbers();
	std::vector<Vector>   pos=pdb.getPositions();
	unsigned k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdb.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdb.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdb.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<pos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// dump the position
	myfile.open ("position.pdb");		
	at=pdbrun.getAtomNumbers();
	std::vector<Vector>   runpos=pdbrun.getPositions();
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<runpos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// now do the alignment
	Tensor Rotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
	std::vector<Vector> alignedpos;
	std::vector<Vector> centeredpos;
	std::vector<Vector> centeredref;
	std::vector<Vector> ddistdpos;
	rmsd->calc_PCA( run, squared, Rotation , ddistdpos, DRotDPos , alignedpos ,centeredpos, centeredref ); 
	myfile.open ("aligned.pdb");		
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<alignedpos[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
	// dump the aligned	
	myfile.open ("reference_centered.pdb");		
	k=0;
	for(std::vector<AtomNumber>::iterator i=at.begin(); i!=at.end(); i++){
		myfile<<"ATOM";
                myfile.width(7);myfile<<std::right<<(*i).serial()<<" "; 
                myfile.width(4);myfile<<std::left<<pdbrun.getAtomName(*i); 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueName(*i)<<" A"; 
                myfile.width(4);myfile<<std::right<<pdbrun.getResidueNumber(*i)<<"    "; 
		myfile.setf( std::ios::fixed, std:: ios::floatfield );
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][0]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][1]*10; 
                myfile.width(8);myfile.precision(3); myfile<<std::right<<centeredref[k][2]*10<<"  1.00  1.00\n"; 
		k++;	
	}	
	myfile.close();			
  }
  // Task 10: derivative of the rotation matrix (in case of reverse transition) 
  if(std::find(task.begin(), task.end(), 10)!=task.end()){
	cout<<"Task 5: calculates the finite difference for derivative of the rotation matrix respect to the the running frame"<<endl;
	Tensor Rotation,OldRotation;
	Matrix<std::vector<Vector> > DRotDPos(3,3);
        std::vector<Vector> DDistDRef;
        std::vector<Vector> alignedpos;
        std::vector<Vector> centeredpos;
        std::vector<Vector> centeredref;
	std::vector<Vector> ddistdpos;
	rmsd->calc_PCA( run, squared, OldRotation , ddistdpos, DRotDPos , alignedpos ,centeredpos, centeredref ); 
	std::vector<Vector> run_save=run;	
	for(unsigned int a=0;a<3;a++){	
		for(unsigned int b=0;b<3;b++){	
			for(unsigned int comp=0;comp<3;comp++){	
				for(unsigned int i=0;i<run.size();++i){
					//change the position		
					run[i][comp]+=eps;
					rmsd->calc_PCA( run, squared, Rotation , ddistdpos, DRotDPos , alignedpos ,centeredpos, centeredref ); 
					cout<<"DROT_DPOS_INVERSE_TRANSFORM COMPONENT "<<comp<<" "<<(Rotation[a][b]-OldRotation[a][b])/(run[i][comp]-run_save[i][comp])<<" "<<DRotDPos[a][b][i][comp]<<"\n";
					// restore the old position
					run=run_save;
				}
			}
		}
	}
  }
	
  return 0;
}