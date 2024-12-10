#include<iostream>
#include<cmath>
#include<numbers>
#include<iomanip>
using namespace std;
// 1. finding a root to x^2-8x+11=0, where x0=6

double f(double x){
    return (x*x - 8*x + 11);
}

double f1(double x){ //first derivative of f
    return x+x-8;
}
// 2. finding length of chord AB
// if 2*theta is the angle AOB, 2*sin(theta)*2 = (2pi-2*theta) is our equation
// 2sin(x) - pi + x = 0
// checking convergence condition for this function gives -0.7<x<1.6
double g(double x){
    return 2*sin(x) + x - M_PI;
}

double g1(double x){ //first derivative of g
    return 2*cos(x) + 1;
}
int count = 0;

double guess(float x0, double err, double (*f)(double x), double(*f1)(double x)){
    if (count > 4000){ // to avoid infinite loop due to too small error
        cout<<"Iteration limit exceeded"<<endl;
        return x0;
    }
    count++;
    double new_guess = x0 - (f(x0) / f1(x0));
    if (abs(new_guess - x0) < err){
        return new_guess;}
    return guess(new_guess, err, f, f1);

}
int main(){
    double x0 = 6;
    double err = pow(10, -5);
    double root = guess(x0, err, &f, &f1);
    cout<<"A root of the equation x^2-8x+11=0 is "<<fixed<<setprecision(5)<<root<<endl;
    x0 = 1; // guessing x0 to be one for ex2
    double err1 = 1e-18; //
    double ans = guess(x0, err1, &g, &g1); // this is the angle theta, lies in convergence area
    cout<<"Length of chord AB is "<<fixed<<setprecision(18)<<sin(ans) * 2<<endl;
}
