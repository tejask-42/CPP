#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;

float f(float x){
    return exp(2*x) + x - 5;
}

int iter = 0;
float guess(float x0, float x1, float(*f)(float x), double err = 1e-5){
    if (iter > 4000){
        cout << "Iteration limit exceeded" << endl;
        return x1;
    }
    iter++;
    float x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));    
    if (abs(x2 - x1) < err){
        return x2;
    }
    return guess(x1, x2, f, err);
}

int main(){
    float x0 = 0, x1 = 1; // f(0) < 0, f(2) > 0
    float ans = guess(x0, x1, &f);
    cout << "A root of e^(2x) + x - 5 is " << fixed << setprecision(5) << ans << endl;
}
