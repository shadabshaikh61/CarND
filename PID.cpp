#include "PID.h"

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp_, double Ki_, double Kd_) {
	Kp = Kp_;
	Ki = Ki_;
	Kd = Kd_;
	
	p_error = 0;
	i_error = 0;
	d_error = 0;
	
}

void PID::UpdateError(double cte) {
	//p_error = Kp*cte;
	//i_error = i_error + cte*Ki;
	//d_error = d_error - cte*Kp;
	p_error = cte;
	i_error = i_error + cte;
	d_error = d_error - cte;
	
}

double PID::TotalError() {
	//return -1*(d_error+i_error+p_error);
	return -1*(Kd*d_error+Ki*i_error+Kp*p_error)/10;
}

