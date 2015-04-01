/***
	Problem: MD simulation of Argon gas, It has the features, (i) periodic boundary
	conditions, (ii) zero center of mass momentum, (iii) a method to fix the temperature 
	Author: T Mathialakan,	A48554709,	MSU.
	Date: March 12, 2015
***/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream> 
#include <string.h> 

#define ATOMS 4

int N = 64;  //Atoms
double L = pow(N, 1.0/3); // N = L^3
double T = 1.0;	// Temperature
double *x = (double *) malloc(sizeof(double)*3*N);	// Position 
double *v = (double *) malloc(sizeof(double)*3*N);	// Velocity 
double *a = (double *) malloc(sizeof(double)*3*N);	// Acceleration
double PI=3.141592653589793;
double pdf =0.0;

void initial_structure();	
void initial_velocity();	
void rescale_velocity();
void print_vector(double *a, int m, int n);

/*
	Initialize FCC lattice 
*/
void initial_structure() { 
	int Li = 1; 
	while (4*pow(Li, 3)< N)	Li++;
	//printf("Li : %d	\n", Li);
	double h = L / Li;	
	// 4 atomic positions in FCC unit cell 
	double x_point[4] = {0.0, 0.5, 0.5, 0.0}; // x coordinate
	double y_point[4] = {0.0, 0.5, 0.0, 0.5}; // y coordinate
	double z_point[4] = {0.0, 0.0, 0.5, 0.5}; // z coordinate
	int Ni = 0;	
	for (int rx = 0; rx < Li; rx++) 
		for (int ry = 0; ry < Li; ry++) 
			for (int rz = 0; rz < Li; rz++) 
				for (int k = 0; k < 4; k++) 
					if (Ni < N) { 
						x[Ni*3] = (rx + x_point[k]) * h; 	
						x[Ni*3+1] = (ry + y_point[k]) * h; 
						x[Ni*3+2] = (rz + z_point[k]) * h;
						Ni++;	
					}
	initial_velocity();
}

/*
	Get random number
*/
double my_rand(){
	double v = 2.0 * rand() / double(RAND_MAX) - 1.0; 
}

/*
	Get random number with a Gaussian probability distribution
*/
double rand_gaussian () { 
	static bool exist = false; 
	static double grand; 
	double f, vs, v1, v2; 
	if (!exist) { 
		do { 
			v1 = my_rand(); 
			v2 = my_rand(); 
			vs = v1 * v1 + v2 * v2; 
			} while (vs >= 1.0 || vs == 0.0); 
				f = sqrt(-2.0 * log(vs) / vs);  
				exist = true; 
				return v2*f;
	} else { 
		exist = false; 
		return v1 * f;
	} 
}

/*
	Calculate the centre of mass velocity in order to make zero centre of mass momentum
*/
double* center_of_mass_velocity(){
	double* vcm = (double *) malloc(sizeof(double)*3);
	for (int j = 0; j < 3; j++)
		vcm[0] = 0; 
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < 3; j++)
			vcm[j] += v[i*3+j];
	for (int j = 0; j < 3; j++)
		vcm[j] /= N; 	
	return vcm;
}

/*
	Compute the initial velocity using Maxwell-Boltzmann distribution
*/
void initial_velocity() {
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < 3; j++) 
			v[i*3+j] = rand_gaussian();
	double* vcm = center_of_mass_velocity();
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < 3; j++){
			v[i*3+j] -= vcm[j];
			}
}

/*
	Calculate the scaling factor lambda 
*/
double cal_lambda(){
	double vs = 0; 
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < 3; j++) 
			vs += v[i*3+j] * v[i*3+j];
	return sqrt( 3*(N-1)*T/vs );
}

/*
	Scale the velocity by the scaling factor lambda
*/
void rescale_velocity() { 	
	double lambda = cal_lambda(); 
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < 3; j++) 
			v[i*3+j] *= lambda;
}

/*
	Calculate the acceleration
*/
void cal_acceleration() {
	pdf =0;
	for (int i = 0; i < N; i++)	
		for (int k = 0; k < 3; k++) 
			a[i*3+k] = 0;
	// Considering all adjacent atoms using periodic boundary condition
	for (int i = 0; i < N-1; i++)	
		for (int j = i+1; j < N; j++) { 
			double xadj[3];	
			double xs = 0;
			for (int k = 0; k < 3; k++) {
				xadj[k] = x[i*3+k] - x[j*3+k];
				if (abs(xadj[k]) > 0.5*L) { 
					if (xadj[k] > 0) 
						xadj[k] -= L;
				else xadj[k] += L;
				}
				xs += xadj[k] * xadj[k];
			}
			double f = 24*(2*pow(xs, -7) - pow(xs, -4)); 
			//pdf+=exp(pow(xs, 2)/-2)/pow(2*PI, 1.5);
			pdf+=xs/(4*N*PI*pow(xs, 2));
			for (int k = 0; k < 3; k++) { 
				a[i*3+k] += xadj[k]*f; 
				a[j*3+k] -= xadj[k]*f;
			}
		}
}

/*
	Calculate velocity using verlet integration algorithm
*/
void verlet_velocity(double dt) {
	cal_acceleration(); // calculate ai(t)
	for (int i = 0; i < N; i++) 
		for (int k = 0; k < 3; k++) { 
			x[i*3+k] += v[i*3+k]*dt + 0.5*a[i*3+k]*dt*dt;
		// Use periodic boundary conditions 
		if (x[i*3+k] < 0) x[i*3+k] += L; 
		if (x[i*3+k] >= L) x[i*3+k] -= L;
		v[i*3+k] += 0.5*a[i*3+k]*dt; // add 0.5 ai(t)dt
	} 
	cal_acceleration();  // calculate ai(t+dt)
	for (int i = 0; i < N; i++) 
		for (int k = 0; k < 3; k++) 
			v[i*3+k] += 0.5*a[i*3+k]*dt; // add 0.5 ai(t+dt)dt
}

/*
	Calculate energy
*/
double cal_energy() { 
	double sum = 0; 
	for (int i = 0; i < N; i++) 
		for (int k = 0; k < 3; k++) 
			sum += v[i*3+k] * v[i*3+k];
	return sum/2;
}

/*
	Calculate temperature
*/
double cal_temperature() { 
	double sum = 0; 
	for (int i = 0; i < N; i++) 
		for (int k = 0; k < 3; k++) 
			sum += v[i*3+k] * v[i*3+k];
	return sum/(3*(N-1));
}

void print_vector(double *a, int m, int n){ 
	int i,j; 
	for(i = 0; i< m; i++){
		for(j = 0; j< n; j++)
			printf("%0.6f\t", a[m*i+j]);
		printf("\n");
	}
	printf("\n");
}     

int main() { 
	initial_structure();
	 double dt = 0.01;  
	 printf("# Temperature \t Temp_diff \t Eng_diff \t Heat_cap \t\n");
	 verlet_velocity(dt); 
	 double T0 = cal_temperature();
	 double E0 = cal_energy();
	 double T1, E1, C;
	 for (int i = 0; i < 1000; i++) {
		T1 = cal_temperature();
		E1 = cal_energy();
		C = (E1-E0)/(2*(T1-T0));
		printf("%f \t %f \t %f \t %f \t %f \n", T0, (T1-T0) , (E1-E0), C, pdf);
		if (i % 200 == 0)
			rescale_velocity();
		verlet_velocity(dt); 
		T0 = T1;
		E0 = E1;
	 }
}