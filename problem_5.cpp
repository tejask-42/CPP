#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdio>
using namespace std;

float f(float x, float y){ // y' = f(x, y)
    return x*x*x*exp(-2*x) - 2*y;
}
float soln(float x){
    return ((pow(x, 4) / 4) + 1) * exp(-2*x);
}

int main(){
    FILE* gp = popen("gnuplot -persist", "w");
    if (!gp){
        cerr<<"Error opening pipe to Gnuplot"<<endl;
        return 1;
    }
    fprintf(gp, "set title 'Runge-Kutta method'\n");
    fprintf(gp, "set xlabel 'X-axis'\n");
    fprintf(gp, "set ylabel 'Y-axis'\n");
    fprintf(gp, "plot '-' using 1:2 with lines title 'step=0.1', '-' using 1:2 with lines title 'Actual function'\n");
    float x0 = 0;
    float y0 = 1;
    float h = 0.1;
    for (int i=0; i<11; i++){
        fprintf(gp, "%f %f\n", x0, y0);
        cout<<"x: "<<x0<<" guess: "<<y0<<" actual: "<<soln(x0)<<endl;
        float k1 = f(x0, y0);
        float k2 = f(x0 + (h/2), y0 + (h*k1/2));
        float k3 = f(x0 + (h/2), y0 + (h*k2/2));
        float k4 = f(x0 + h, y0 + h*k3);
        y0 = y0 + ((h / 6) * (k1 + 2*k2 + 2*k3 + k4));
        x0 = x0 + h;
    }
    fprintf(gp, "e\n");
    for (int k=0; k<1001; k++){
        float x = float(k) / 1000;
        float y = soln(x);
        fprintf(gp, "%f %f\n", x, y);
    }
    fprintf(gp, "e\n");
    pclose(gp);
}