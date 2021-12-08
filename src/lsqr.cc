#include <limits>
#include <iostream>

#include "lsqr.h"

/********************************************************************************
*Few of the code blocks mentioned below are taken from open source git repo
* https://github.com/harusametime/LSQRwithEigen/blob/master/LSQR.cpp
*********************************************************************************/

LSQR::LSQR(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const double eps)
{
    mat_a = A;
    vec_b = b;
    this->eps = eps;

    assert( mat_a.rows() > 1 );
    assert( mat_a.cols() > 1 );

    assert( mat_a.rows() == vec_b.size() );
}

Eigen::MatrixXd LSQR::SolveForX()
{
	/***********************
	* Variables initialisation
	************************/
	vec_x.setZero(mat_a.cols(), 1);
	double beta = vec_b.norm();
	Eigen::VectorXd vec_u = vec_b / beta;
	Eigen::VectorXd temp_vec = mat_a.transpose() * vec_u;
	double alpha = temp_vec.norm();
	Eigen::VectorXd vec_v = temp_vec/alpha;
	Eigen::VectorXd vec_w = vec_v;
	double phi_bar = beta;
	double rho_bar = alpha;
	int max_itrs = 1000;

	/*************************************************************
	* Initialisation of variables used for the test of convergence
	**************************************************************/
	double z = 0;
	double cs2 = -1;
	double sn2 = 0;
	double dd_norm = 0;
	double b_norm = beta;
	double r_norm = beta;
	double x_norm = 0;
	double xx_norm = 0;
	double mat_a_norm = 0;
	double mat_a_cond = 0;


    int itr = 0;
	while (itr < max_itrs){
		
		/*************************************
		* Bidiagnolization of matrix
		**************************************/
		Eigen::VectorXd rhs_beta = mat_a *vec_v - alpha * vec_u;
		beta = rhs_beta.norm();
		vec_u = rhs_beta / beta;

		Eigen::VectorXd rhs_alpha = mat_a.transpose() * vec_u  - beta * vec_v;
		alpha = rhs_alpha.norm();
		vec_v = rhs_alpha / alpha;

		/***************************************************************************
		* Construct and apply next orthogonal transformation to bidigonalized matrix
		****************************************************************************/

		double rho = sqrt(rho_bar * rho_bar + beta * beta);
		double c = rho_bar / rho;
		double s = beta / rho;
		double theta = s * alpha;
		rho_bar = -c* alpha;
		double phi = c * phi_bar;
		phi_bar = s*phi_bar;


		/****************************************************
		* Tests to check if the algorithm reached convergence
		*****************************************************/

		double gamma_bar = -cs2 *rho;
		double rhs = phi - sn2 * rho * z;
		double zbar = rhs / gamma_bar;
		x_norm = sqrt(xx_norm + zbar * zbar);
		double gamma = sqrt(gamma_bar* gamma_bar + theta* theta);
		cs2 = gamma_bar / gamma;
		sn2 = theta / gamma;
		z = rhs / gamma;
		xx_norm += z * z;


		Eigen::VectorXd rhow = (1 / rho) * vec_w;
		dd_norm = dd_norm + rhow.norm() * rhow.norm();
		mat_a_norm = sqrt(mat_a_norm * mat_a_norm + alpha * alpha + beta * beta);
		mat_a_cond = mat_a_norm + sqrt(dd_norm);
		r_norm = phi_bar;

		double Arnorm = alpha * abs(s * phi);
		double test1 = r_norm / b_norm;
		double test2 = 0;
		double test3 = 0;

		if (mat_a_norm == 0 || r_norm == 0){
			test2 = std::numeric_limits<double>::max();
		}
		else{
			test2 = Arnorm / (mat_a_norm * r_norm);
		}
		if (mat_a_cond == 0){
			test3 = std::numeric_limits<double>::max();
		}
		else{
			test3 = 1 / mat_a_cond;
		}
		double t1 = test1 / (1 + mat_a_norm*x_norm / b_norm);
		double rtol = eps + eps * mat_a_norm * x_norm / b_norm;

		
		
		if (test1 <= rtol || test2 <= eps || test3 <= eps){
			break;
		}


		/*************************************
		* Updates the vectors x and w
		**************************************/
		vec_x = vec_x + (phi / rho) * vec_w;
		vec_w = vec_v - (theta / rho) * vec_w;

		itr++;
    }

    return vec_x;
}
