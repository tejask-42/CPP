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
    fprintf(gp, "set title 'Euler''s method'\n");
    fprintf(gp, "set xlabel 'X-axis'\n");
    fprintf(gp, "set ylabel 'Y-axis'\n");
    fprintf(gp, "plot '-' using 1:2 with lines title 'h=0.1', '-' using 1:2 with lines title 'h=0.05', '-' using 1:2 with lines title 'h=0.025', '-' using 1:2 with lines title 'Actual function'\n");

    float steps[3] = {0.1, 0.05, 0.025};
    float x_vals[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for (int i=0; i<3; i++){
        float step = steps[i];
        int runs = 1 / step;
        cout<<"Applying Euler's method with step size "<<step<<":"<<endl;
        float x0 = 0;
        float y0 = 1; // initial condition
        for (int j=0; j<runs; j++){
            float m = f(x0, y0);
            fprintf(gp, "%f %f\n", x0, y0);
            float y1 = y0 + step * m;
            float x1 = x0 + step;
            for (int k=0; k < 10; k++){ // checking values for 0.1, 0.2 etc
                if (abs(x_vals[k] - x1) < 1e-2){
                    cout<<"x: "<<x1<<" guess: "<<y1<<" actual: "<<soln(x1)<<endl;
                }
            }
            x0 = x1;
            y0 = y1;
        }
        fprintf(gp, "e\n");
        cout<<"--------------------------------------------------"<<endl;
    }
    for (int k=0; k<1001; k++){
        float x = float(k) / 1000;
        float y = soln(x);
        fprintf(gp, "%f %f\n", x, y);
    }
    fprintf(gp, "e\n");
    pclose(gp);
}
