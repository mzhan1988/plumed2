#include "HistogramBead.h"
#include <vector>
#include <limits>
#include "Tools.h"

using namespace PLMD;

std::string HistogramBead::documentation( bool dir ) {
  std::ostringstream ostr;
  if(dir){
     ostr<<"\\f$ w(d)=\\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma_d} \\exp\\left( -\\frac{d^2}{2\\sigma_d^2} \\right) \\textrm{d}d \\f$";
     ostr<<"where \\f$ \\sigma_d=(b_d-a_d)k_d \\f$.  The parameters of the functions are specifed in fractional coordinates using ";
     ostr<<"(XLOWER=\\f$a_x\\f$ XUPPER=\\f$b_x\\f$ XSMEAR=\\f$k_x\\f$ YLOWER=\\f$a_y\\f$ YUPPER=\\f$b_y\\f$ YSMEAR=\\f$k_y\\f$ ZLOWER=\\f$a_z\\f$ ZUPPER=\\f$b_z\\f$ ZSMEAR=\\f$k_z\\f$).";
     ostr<<"If any of the SMEAR keywords are not present then the default \\f$k_x=0.5\\f$ is used in that direction. ";
  } else {
     ostr<<"\\f$ w(r)=\\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma} \\exp\\left( -\\frac{r^2}{2\\sigma^2} \\right) \\textrm{d}r \\f$";
     ostr<<"where \\f$ \\sigma=(b-a)k \\f$.  The parameters of the function are specified using (LOWER=\\f$a\\f$ UPPER=\\f$b\\f$ SMEAR=\\f$k\\f$). ";
     ostr<<"If the SMEAR keyword is not present then by default \\f$k=0.5\\f$.";
  }
  return ostr.str();
}

std::string HistogramBead::description() const {
  std::ostringstream ostr;
  ostr<<"betweeen "<<lowb<<" and "<<highb<<" width of gaussian window equals "<<width;
  return ostr.str();
}

std::string HistogramBead::histodocs( bool dir ) {
  std::ostringstream ostr;
  if(dir){
     ostr<<"The range is divided into a discete number of bins and the number of values that fall within each bin is calculated using ";
     ostr<<"\\f$ w(x)=\\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma_x} \\exp\\left( -\\frac{x^2}{2\\sigma_x^2} \\right) \\textrm{d}x \\f$";
     ostr<<"where \\f$ \\sigma_x=(b_x-a_x)*k_x \\f$.  The particular range of interest and number of bins are specified using ";
     ostr<<"(XNBINS=\\f$n\\f$ XLOWER=\\f$a_x\\f$ XUPPER=\\f$b_x\\f$ XSMEAR=\\f$k_x\\f$). If the SMEAR keyword is not present then by default \\f$k=0.5\\f$";
  } else { 
     ostr<<"The range is divided into a discete number of bins and the number of values that fall within each bin is calculated using ";
     ostr<<"\\f$ w(r)=\\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma} \\exp\\left( -\\frac{r^2}{2\\sigma^2} \\right) \\textrm{d}r \\f$";
     ostr<<"where \\f$ \\sigma=(b-a)k \\f$.  The particular range of interest and number of bins are specified using ";
     ostr<<"(NBINS=\\f$n\\f$ LOWER=\\f$a\\f$ UPPER=\\f$b\\f$ SMEAR=\\f$x\\f$). If the SMEAR keyword is not present then by default \\f$k=0.5\\f$";
  }
  return ostr.str();
}

void HistogramBead::generateBins( const std::string& params, const std::string& dd, std::vector<std::string>& bins ){
  if( dd.size()!=0 && params.find(dd)==std::string::npos) return;
  std::vector<std::string> data=Tools::getWords(params);
  plumed_assert(data.size()>=1);

  unsigned nbins; std::vector<double> range(2); std::string smear;
  bool found_nb=Tools::parse(data,dd+"NBINS",nbins);
  plumed_massert(found_nb,"Number of bins in histogram not found");
  bool found_r=Tools::parse(data,dd+"LOWER",range[0]);
  plumed_massert(found_r,"Lower bound for histogram not specified");
  found_r=Tools::parse(data,dd+"UPPER",range[1]);
  plumed_massert(found_r,"Upper bound for histogram not specified");
  plumed_massert(range[0]<range[1],"Range specification is dubious"); 
  bool found_b=Tools::parse(data,dd+"SMEAR",smear);
  if(!found_b){ Tools::convert(0.5,smear); }  

  std::string lb,ub; double delr = ( range[1]-range[0] ) / static_cast<double>( nbins );
  for(unsigned i=0;i<nbins;++i){
     Tools::convert( range[0]+i*delr, lb );
     Tools::convert( range[0]+(i+1)*delr, ub );
     bins.push_back( dd + "LOWER=" + lb + " " + dd + "UPPER=" + ub + " " + dd + "SMEAR=" + smear );
  }
  plumed_assert(bins.size()==nbins);
  if( dd.size()==0 ) plumed_massert(data.size()==0,"Error reading histogram"); 
}

void HistogramBead::set( const std::string& params, const std::string& dd, std::string& errormsg ){
  if( dd.size()!=0 && params.find(dd)==std::string::npos) return;
  std::vector<std::string> data=Tools::getWords(params);
  plumed_assert(data.size()>=1);

  double smear;
  bool found_r=Tools::parse(data,dd+"LOWER",lowb);
  if( !found_r ) errormsg="Lower bound has not been specified use LOWER";
  found_r=Tools::parse(data,dd+"UPPER",highb);
  if( !found_r ) errormsg="Upper bound has not been specified use UPPER"; 
  if( lowb>=highb ) errormsg="Lower bound is higher than upper bound"; 
  
  smear=0.5; bool found_b=Tools::parse(data,dd+"SMEAR",smear);
  width=smear*(highb-lowb); init=true;
  if( dd.size()==0 ){ if( data.size()!=0) errormsg="Error reading within"; }
}

void HistogramBead::printKeywords( Log& log ) const {
  Keywords hkeys;
  hkeys.add("compulsory","LOWER","the lower boundary for this particular bin");
  hkeys.add("compulsory","UPPER","the upper boundary for this particular bin");
  hkeys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution"); 
  hkeys.print( log );
}