#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	cout << "init\n";
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);
	
	P_ <<   10, 0, 0, 0, 0,
			0, 10, 0, 0, 0,
			0, 0, 1000, 0, 0,
			0, 0, 0, 10, 0,
			0, 0, 0, 0, .01;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.65;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.8;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	 TODO:

	 Complete the initialization. See ukf.h for other member properties.

	 Hint: one or more values initialized above might be wildly off...
	 */

	// time when the state is true, in us
	time_us_ = 0;

	Xsig_pred_ = MatrixXd(5, 15);

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	///* Sigma point spreading parameter
	lambda_ = 3 - n_aug_;
	
	///* Weights of sigma points
	weights_ = VectorXd(15);
	double w0 = lambda_ / (lambda_ + n_aug_);
	double w1 = w0 / (2 * lambda_);
	int n_a = 1 + 2 * n_aug_;
	  
	weights_(0) = w0;
	int i = 1;
	while (i < n_a){
	    weights_(i) = w1;
	    i++;
	}

	is_initialized_ = false;
}

UKF::~UKF() {
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	 TODO:

	 Complete this function! Make sure you switch between lidar and radar
	 measurements.
	 */
	if (!is_initialized_){
		cout << "UKF\n";
		
		if (meas_package.sensor_type_ = MeasurementPackage::RADAR){
		  //cout << "!init radar \n";
		  float ro = meas_package.raw_measurements_(0);
		  float theta = meas_package.raw_measurements_(1);
		  x_(0) = ro * cos(theta);
		  x_(1) = ro * sin(theta);
		  //cout << "x_ \n" << x_ << "\n";
		}
	    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	      /**
	      Initialize state.
	      */
	      //cout << "!init lidar \n";
	      x_(0) = meas_package.raw_measurements_(0);
	      x_(1) = meas_package.raw_measurements_(1);
	      //cout << "x_ \n" << x_ << "\n";
	    }
	time_us_ = meas_package.timestamp_;
	is_initialized_ = true;
	return;
	}
	
	double delta_t;
	delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
	time_us_= meas_package.timestamp_;
	Prediction(delta_t);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
		UpdateRadar(meas_package);
	}
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    	UpdateLidar(meas_package);
    }
	//cout << "\n";
	return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	 TODO:

	 Complete this function! Estimate the object's location. Modify the state
	 vector, x_. Predict sigma points, the state, and the state covariance matrix.
	 */
	 //  Create sigma poinnts 
	 // refer to Term 2, Lesson 7, Section 17.
	 //create augmented mean vector
	 VectorXd x_aug = VectorXd(7);
	 x_aug.fill(0.0);

	 //create augmented state covariance
	 MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	 P_aug.fill(0.0);
	 
	 MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	 x_aug << x_, 0, 0;
	 
	 //create augmented covariance matrix
	 MatrixXd Q = MatrixXd(2,2);
	 Q << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
	 P_aug.block<5, 5>(0,0) = P_;
	 P_aug.block<2, 2>(n_x_, n_x_) = Q;
	 
	  
	 //create square root matrix
	 // MatrixXd A = MatrixXd(n_aug_, n_aug_);
	 // A = P_aug.llt().matrixL();
	 
	 MatrixXd A = MatrixXd(n_aug_, n_aug_);
	 A.fill(0.0);
	 
	 MatrixXd F = P_aug;
		int j = 0;
		while (j < n_aug_){
			int k = 0;
			double sk = 0;
			while (k < j){
				sk += A(j, k) * A(j, k);
				k++;
			}
			A(j, j) = sqrt((F(j, j) - sk));
			int i = j + 1;
			while (i < n_aug_){
				int l = 0;
				double sl = 0;
				while(l < j){
					sl += A(i, l) * A(j, l);
					l++;
				}
				A(i, j) = (F(i, j) - sl)/A(j, j);
				i++;
			}
			j++;
		}
		
	 //cout << "\n \nPREDICTION\n";
	 //cout << "delta_t: " << delta_t << "\n";
	 //cout << "A \n" << A << "\n";
	 //cout << "A * A' \n" << A * A.transpose() << "\n";
	 //cout << "P_aug \n" << P_aug << "\n"; 

	 //create augmented sigma points
	 Xsig_aug.col(0) = x_aug;

	 MatrixXd dx = MatrixXd(n_aug_, n_aug_);
	 dx = (sqrt(lambda_ + n_aug_)) * A;
	 
	 //cout << "lambda _ n_aug \n" << lambda_ + n_aug_ << "\n";
	 //cout << "dx \n" << dx << "\n";
	 
	 int i = 1;
	 while(i<=n_aug_){
	     Xsig_aug.col(i)=x_aug + dx.col(i-1);
	     i++;
	 }
	   while(i<=n_aug_*2){
	     Xsig_aug.col(i) = x_aug - dx.col(i-1-n_aug_);
	     i++;
	 }
	 
	 // predict sigma points
	 // refer to term 2, lesson 7, section 20
	 
	 //cout << "Xsig_aug \n" << Xsig_aug << "\n";
	 
	 Xsig_pred_.fill(0.0);
	 i = 0;
	 while (i <n_aug_*2 + 1){
	     VectorXd x = VectorXd(n_aug_);
	     VectorXd xo = VectorXd(n_x_);
	     VectorXd dx1 = VectorXd(n_x_);
	     VectorXd dx2 = VectorXd(n_x_);
	     x = Xsig_aug.col(i);
	     if (fabs(x(4)) > 0.0001){
	         dx1 << x(2)/x(4)*( sin(x(3)+x(4)*delta_t)-sin(x(3))),
	                x(2)/x(4)*(-cos(x(3)+x(4)*delta_t)+cos(x(3))),
	                0,
	                x(4)*delta_t,
	                0;
	         //cout << "thd = 0 dx1 \n" << dx1 << "\n";
	         dx2 << 0.5*delta_t*delta_t*cos(x(3))*x(5),
	                0.5*delta_t*delta_t*sin(x(3))*x(5),
	                delta_t*x(5),
	                0.5*delta_t*delta_t*x(6),
	                delta_t*x(6);
	         //cout << "thd = 0 dx2 \n" << dx2 << "\n";
	         xo = x.head(n_x_) + dx1 + dx2;
	     }
	     else{
	         dx1 << x(2)*cos(x(3))*delta_t,
	                x(2)*sin(x(3))*delta_t,
	                0,
	                x(4)*delta_t,
	                0;
	         //cout << "thd <> 0 dx1 \n" << dx1 << "\n";
	         dx2 << 0.5*delta_t*delta_t*cos(x(3))*x(5),
	                0.5*delta_t*delta_t*sin(x(3))*x(5),
	                delta_t*x(5),
	                0.5*delta_t*delta_t*x(6),
	                delta_t*x(6);
	         //cout << "thd <> 0 dx1 \n" << dx2 << "\n";
	         xo = x.head(n_x_) + dx1 + dx2;
	     }
	     Xsig_pred_.col(i) = xo;
      i++;
	  }
	 
	  // predict xbar, cov(x)
	  // refer to term 2, lesson 7, section 23

	  //predict state mean
	  i = 0;
	  x_.fill(0.0);
	  int n_a = 1 + 2 * n_aug_;
	  while (i < n_a){
	      x_ += weights_(i) * Xsig_pred_.col(i);
	      i++;
	  }
	  
	  //predict state covariance matrix
	  i = 0;
	  P_.fill(0.0);
	  while (i < n_a){
	    // state difference
	    VectorXd x_diff = Xsig_pred_.col(i) - x_;
	    //angle normalization
	    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	    P_ += weights_(i) * x_diff * x_diff.transpose();
	    i++;
	  }

		 //cout << "Xsig_pred_ \n" << Xsig_pred_ << "\n";
		 //cout << "x_ \n" << x_ << "\n";
		 //cout << "P_ \n" << P_ << "\n";
	 return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	 TODO:

	 Complete this function! Use lidar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.

	 You'll also need to calculate the lidar NIS.
	 */
	MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
	MatrixXd H = MatrixXd(2, n_x_);
	H << 1, 0, 0, 0, 0,
		 0, 1, 0, 0, 0;
	
	MatrixXd R = MatrixXd(2, 2);
	R << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;
	
	VectorXd y = meas_package.raw_measurements_ - H * x_;
	MatrixXd Ht = H.transpose();
	MatrixXd P_Ht = P_ * Ht;
	MatrixXd S = H * P_Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd K =  P_Ht * Si;

	//new state
	x_ = x_ + (K * y);
	P_ = (I - K * H) * P_;
	 //cout << "UPDATE LIDAR \n";
	 //cout << "raw meas \n" << meas_package.raw_measurements_ << "\n";
	 //cout << "x_ \n" << x_ << "\n";
	 //cout << "P_ \n" << P_ << "\n\n";
	return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	 TODO:

	 Complete this function! Use radar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.

	 You'll also need to calculate the radar NIS.
	 */
	 // refer to term 2, lesson 7, section 26
	 // transform sigma points into measurement space
	 int i = 0;
	 int n_a = 1 + 2 * n_aug_;
	 //create matrix for sigma points in measurement space
	 MatrixXd Zsig = MatrixXd(3, 15);
	 Zsig.fill(0.0);

	 //mean predicted measurement
	 VectorXd z_pred = VectorXd(3);
	 z_pred.fill(0.0);
	  
	 //measurement covariance matrix S
	 MatrixXd S = MatrixXd(3, 3);
	 S.fill(0.0);

	 while (i < n_a){
	     VectorXd x = Xsig_pred_.col(i);
	      
	     double px = x(0);
	     double py = x(1);
         double v = x(2);	      
         double psi = x(3);
	     double rho = 0;
	     double theta = 0;
	     double ro_dot = 0;
	     if (fabs(px) < 0.0001 && fabs(py) < 0.0001){
	       px = 0.0001;
	       py = 0.0001;
	     }
	     
	     rho = sqrt(px*px + py*py);
	     if (fabs(px) < 0.0001){
	    	 px = 0.0001;
	     }
	     theta = atan2(py, px);
	     ro_dot = (px *  cos(psi) * v + py * sin(psi) * v) / rho;
	      
	     Zsig.col(i) << rho, theta, ro_dot;
	     i++;
	 }
	 //calculate mean predicted measurement
	 i = 0;
	 while(i < n_a){
	     z_pred += weights_(i) * Zsig.col(i); 
	     i++;
	 }
	  
	 //calculate measurement covariance matrix S
	 MatrixXd R = MatrixXd(3, 3);
	 R.fill(0.0);
	 R(0,0) = std_radr_*std_radr_;
	 R(1,1) = std_radphi_*std_radphi_;
	 R(2,2) = std_radrd_*std_radrd_;

	 //cout << "R \n" << R << "\n";
	 i = 0;
	 VectorXd z_diff = VectorXd(3);
	 z_diff.fill(0.0);
	 while(i < n_a){
	     z_diff = Zsig.col(i) - z_pred;
	     //angle normalization
	     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
	     S += weights_(i) * z_diff * z_diff.transpose();
	     i++;
	 }
	 S += R;

	 //cout << "S \n" << S << "\n";


	 // refer to term 2, lessson 7, section 29
	 //calculate cross correlation matrix

	 MatrixXd T = MatrixXd(5, 3);
	 T.fill(0.0);

	 i = 0;
	 while (i < n_a){
	   VectorXd x_diff = VectorXd(5);
	   x_diff = Xsig_pred_.col(i) - x_;
	   while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	   while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	   T = T + weights_(i) * x_diff * z_diff.transpose();
	   i++;
     }

	 //cout << "T \n" << T << "\n";
	 //cout << "Sinv \n" << S.inverse() << "\n";
	 //calculate Kalman gain K;
	 MatrixXd K = MatrixXd(5, 3);
	 K = T * S.inverse();
	  
	 //cout << "K \n" << K << "\n";
	 //update state mean and covariance matrix
	 VectorXd zpred_diff = VectorXd(3);
	 
	 VectorXd z = meas_package.raw_measurements_;

	 //cout << "sensor type\n" << meas_package.sensor_type_ << "\n";
	 //cout << "z\n" << z << "\n";
	 //cout << "z_pred\n" << z_pred << "\n";
	 
	 zpred_diff = z - z_pred;
	 //angle normalization
	 while (zpred_diff(1)> M_PI) zpred_diff(1)-=2.*M_PI;
	 while (zpred_diff(1)<-M_PI) zpred_diff(1)+=2.*M_PI;
	    
	 x_ = x_ + K * zpred_diff;
	 P_ = P_ - K * S * K.transpose();
	 //cout << "UPDATE RADAR \n";
	 //cout << "raw meas \n" << meas_package.raw_measurements_ << "\n";
	 //cout << "x_ \n" << x_ << "\n";
	 //cout << "P_ \n" << P_ << "\n";
	 return;
}
