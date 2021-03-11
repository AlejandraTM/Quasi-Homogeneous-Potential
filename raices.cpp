

/* How to run the code in POP_UP! OS or Ubuntu:
g++ raices.cpp -o raices -lgsl -lgslcblas -lm
./raices
*/

#include <iostream> 
#include <gsl/gsl_roots.h>
#include <fstream>
#include <stdlib.h>

/*Parameters definition and function for relative equilibium*/
struct Params {double C,a,b,alpha,beta;};
double examplefunction(double x, void* param);
/*This function calculates the roots in an interval from x_lo to x_hi using the GSL library*/
double roots(double x_lo, double x_hi, void* param);
/*This function determines the root stability*/
double determinant(double root, void* param, double A1, double B1);

int main() {
/*Constan values declaration*/
 double M=1.0;
 double A=10.0;
 double B=1.0;
 double A1=10.0;
 double B1=1.0;
 double alpha=M*(A+4.0*A1);
 double beta=M*(B+16.0*B1);
 double a=4.0;
 double b=8.0;
/*Function limits*/
 double x_lo;
 double x_hi;
/*Root counter*/
 double num;
/*Root and determinant*/
double r,d; 
 
/*Open Files Code*/
 std::ofstream archivo2;
 archivo2.open("Roots2-(8,4).txt",std::ios::out);
 if(archivo2.fail())
	{
		std::cout<<"No se pudo abrir el archivo."<<"\n";
		exit(1);
	}
	
/*Roots calculation*/
for(double C=0.0;C<=10.0;C=C+0.1)//Values of angular momentum from -1 to 1
{
	//archivo2<<"Potential values: "<<a<<" and "<<b<<" Angular momenta: "<<C<<"\n";
	num=0.0;
	Params args = {C,a,b,alpha,beta};
	//std::cout<<"momentum:"<<C<<"\n";
	for(double j=-100.0;j<100.0;j=j+0.01)//Divition of the function interval
	{
		x_lo=j-0.01;
		x_hi=j;
     		if((examplefunction(x_lo,&args)*examplefunction(x_hi,&args))<0.0)
		{
			num=num+1.0;//Number roots counter
			r=roots(x_lo,x_hi,&args);//Roots calculation using the roots function
			d=determinant(r,&args,A1,B1);//Stability calculation using the Hessian determinant condition
			if(d==0.0){
				archivo2<<C<<";"<<r<<";"<<"Not results"<<"\n";
			}
			else{
				if(d<0.0){
					archivo2<<C<<";"<<r<<";"<<"Unstable"<<"\n";
				}
				else{
					archivo2<<C<<";"<<r<<";"<<"Stable"<<"\n";
				}
			}
		}
	}
	//std::cout<<"numero:"<<num<<"\n";
}

 /*Closed file codes*/
	archivo2.close();
	return 0;
}

/*Functions declarated*/
double examplefunction(double x, void* param){
	Params* p = (Params*)param;
	return -2.0*pow(p->C,2.0)*pow(x,(p->a+p->b))+pow(x,2.0)*(pow(x,p->b)*(p->a)*(p->alpha)+pow(x,p->a)*(p->b)*(p->beta));
}

double roots(double x_lo, double x_hi, void* param){
	Params* p = (Params*)param;
	double raiz=0.0;
	Params args = {p->C,p->a,p->b,p->alpha,p->beta};
	gsl_root_fsolver* solver;
	gsl_function      fwrapper;
	solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	fwrapper.function = examplefunction;
	fwrapper.params = &args;
	gsl_root_fsolver_set(solver,&fwrapper,x_lo,x_hi);
	int status = 1;
	for (int iter=0; status and iter < 10; ++iter){
		gsl_root_fsolver_iterate(solver);
		double x_rt = gsl_root_fsolver_root(solver);
		double x_lo = gsl_root_fsolver_x_lower(solver);
		double x_hi = gsl_root_fsolver_x_upper(solver);
		raiz=x_rt;
		status=gsl_root_test_interval(x_lo,x_hi,0,1e-6);
	}
	gsl_root_fsolver_free(solver);
	return raiz;
}

double determinant(double root, void* param, double A1, double B1){
	Params* p = (Params*)param;
	double determinant;
	double f1,f2;//Those funcions take M=1. It it is different f1 and f2 are different too.
	f1=-((p->a*(p->a+1.0)*p->alpha)/(pow(root,p->a+2.0)))-((p->b*(p->b+1)*p->beta)/(pow(root,p->b+2.0)))+((6.0*pow(p->C,2.0))/(1.0*pow(root,4)));
	f2=((16.0*p->a*A1)/(pow(root,p->a+2.0)))+((64.0*p->b*B1)/(pow(root,p->b+2.0)));
	determinant=f1*f2;
	//std::cout<<f1<<" ; "<<f2<<"\n";
	return determinant;
}
