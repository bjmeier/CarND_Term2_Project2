#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size()==0){
	  cout << "no estimations" << endl;  
	  return rmse;
	}
	if( estimations.size() != ground_truth.size()){
      cout << "inconsistant estimation and ground truth points" << endl;
	  return rmse;
	}
	for(int i=0; i < estimations.size(); ++i){
	  if (estimations[i].rows() != ground_truth[i].rows()){
        cout << "inconsistant estimation and ground truth rows" << endl;
	    return rmse;
	  }
	}
	
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
	  // ... your code here
	  VectorXd res = (estimations[i]-ground_truth[i]).array().square();
	  rmse = rmse + res;
	  
	}
	rmse = rmse.array()/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}