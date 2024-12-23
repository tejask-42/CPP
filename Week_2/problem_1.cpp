#include<iostream>
#include<cmath>
#include<numbers>
#include<iomanip>
#include<cstdio>
using namespace std;

float f(float t, float z1, float z2){ // y''
    return sin(t) + 2*cos(t) - 3*z1 - 4*z2;
}
float soln(float t){ // y
    return 0.5*sin(t) + 0.75*exp(-1*t) + 0.25*exp(-3*t);
}
int main(){
        FILE* gp = popen("gnuplot -persist", "w");
    if (!gp){
        cerr<<"Error opening pipe to Gnuplot"<<endl;
        return 1;
    }
    fprintf(gp, "set title 'Euler''s method and Runge-Kutta method'\n"); 
    fprintf(gp, "set xlabel 'X-axis'\n");
    fprintf(gp, "set ylabel 'Y-axis'\n");
    fprintf(gp, "plot '-' using 1:2 with lines title 'Euler h=0.1', '-' using 1:2 with lines title 'Euler h=0.01', '-' using 1:2 with lines title 'Actual function', '-' using 1:2 with lines title 'RK4 step=0.1'\n");

    float steps[2] = {0.1, 0.01};
    for (int i=0; i<2; i++){
        float step = steps[i];
        int runs = 1 / step;
        float t0 = 0;
        float z1 = 1; // initial condition
        float z2 = -1;
        cout<<"Applying Euler's method with step size "<<step<<":"<<endl;
        for (int j=0; j<runs; j++){
            float m = f(t0, z1, z2);
            fprintf(gp, "%f %f\n", t0, z1);
            float f1 = f(t0, z1, z2);
            t0 = t0 + step;
            z1 = z1 + step*z2;
            z2 = z2 + step*f1;
            cout<<"x: "<<t0<<" guess: "<<z1<<" actual: "<<soln(t0)<<endl;
        }
        cout<<"--------------------------------------------------"<<endl;
        fprintf(gp, "e\n");
    }
    for (int k=0; k<1001; k++){
        float x = float(k) / 1000;
        float y = soln(x);
        fprintf(gp, "%f %f\n", x, y);
    }
    fprintf(gp, "e\n");
    float t0 = 0;
    float y = 1;
    float y1 = -1;
    float h = 0.1;
    cout<<"Applying Runge-Kutta method with step size 0.1:"<<endl;
    for (int i=0; i<11; i++){
        cout<<"x: "<<t0<<" guess: "<<y<<" actual: "<<soln(t0)<<endl;
        fprintf(gp, "%f %f\n", t0, y);
        float k11 = y1;
        float k12 = f(t0, y, y1);
        float k21 = y1 + (h / 2) * k12;
        float k22 = f(t0 + h / 2, y + (h / 2) * k11, y1 + (h / 2) * k12);
        float k31 = y1 + (h / 2) * k22;
        float k32 = f(t0 + h / 2, y + (h / 2) * k21, y1 + (h / 2) * k22);
        float k41 = y1 + h * k32;
        float k42 = f(t0 + h, y + h * k31, y1 + h * k32);

        y += (h / 6) * (k11 + 2 * k21 + 2 * k31 + k41);
        y1 += (h / 6) * (k12 + 2 * k22 + 2 * k32 + k42);
        // Update t
        t0 += h;
    }
    pclose(gp);
}
