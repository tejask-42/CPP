#include <iostream>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <unistd.h>
using namespace std;


float f(float x) {
    return x - (x * x * x) / 3;
}
float x_dot(float x, float y, float mu) {
    return mu * (x - (x * x * x) / 3 - y);
}

float y_dot(float x, float mu) {
    return x / mu;
}

// RK4 Method
void rk4Method(vector<float>& x_rk, vector<float>& y_rk, float h, float mu, float t, float x, float y, int steps) {
    for (int i = 0; i < steps; i++) {
        x_rk.push_back(x);
        y_rk.push_back(y);

        float k1x = x_dot(x, y, mu);
        float k1y = y_dot(x, mu);

        float x_k2 = x + h * k1x / 2;
        float y_k2 = y + h * k1y / 2;
        float k2x = x_dot(x_k2, y_k2, mu);
        float k2y = y_dot(x_k2, mu);

        float x_k3 = x + h * k2x / 2;
        float y_k3 = y + h * k2y / 2;
        float k3x = x_dot(x_k3, y_k3, mu);
        float k3y = y_dot(x_k3, mu);

        float x_k4 = x + h * k3x;
        float y_k4 = y + h * k3y;
        float k4x = x_dot(x_k4, y_k4, mu);
        float k4y = y_dot(x_k4, mu);

        x += h * (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y += h * (k1y + 2 * k2y + 2 * k3y + k4y) / 6;

        t += h;
    }
}

// Euler Method
void eulerMethod(vector<float>& x_eul, vector<float>& y_eul, float h, float mu, float t, float x, float y, int steps) {
    for (int i = 0; i < steps; i++) {
        x_eul.push_back(x);
        y_eul.push_back(y);

        float dx = x_dot(x, y, mu);
        float dy = y_dot(x, mu);
        x += h * dx;
        y += h * dy;
        t += h;
    }
}
float mu = 7.0;



int main() {
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        cerr << "Error opening pipe to Gnuplot" << endl;
        return 1;
    }
    fprintf(gp, "set multiplot layout 2,2 rowsfirst\n");

    float h = 0.1;
    int steps = 200;
    float t = 0.0;
    float x = 1.0;
    float y = 0.0;

    vector<float> x_eul, y_eul, x_rk, y_rk;
    eulerMethod(x_eul, y_eul, h, mu, t, x, y, steps);
    rk4Method(x_rk, y_rk, h, mu, t, x, y, steps);

    vector<float> x_vals; // to use when plotting y = x - x^3/3
    for (int i = -steps; i < steps; i++) {
        x_vals.push_back(0.01 * i);
    }
    float t_ind = 0;
    vector<float> time; // to use for plotting time
    for (int i=0; i<steps; i++){
        time.push_back(t_ind);
        t_ind += h;
    }

    fprintf(gp, "set title 'Van der Pol Oscillator: Euler vs RK4'\n");
    fprintf(gp, "set xlabel 'Time'\n");
    fprintf(gp, "set ylabel 'x'\n");
    fprintf(gp, "plot '-' with lines title 'Euler x', '-' with lines title 'RK4 x'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], x_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], x_rk[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set title 'Van der Pol Oscillator: Euler vs RK4'\n");
    fprintf(gp, "set xlabel 'Time'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "plot '-' with lines title 'Euler y', '-' with lines title 'RK4 y'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], y_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], y_rk[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set title 'y vs x Euler method'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "plot '-' with lines title 'Euler', '-' with lines title 'y=x-x3/3'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", x_eul[i], y_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < 2 * steps; i++) {
        fprintf(gp, "%f %f\n", x_vals[i], f(x_vals[i]));
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set title 'y vs x RK4 method'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "plot '-' with lines title 'RK4', '-' with lines title 'y=x-x3/3'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", x_rk[i], y_rk[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < 2 * steps; i++) {
        fprintf(gp, "%f %f\n", x_vals[i], f(x_vals[i]));
    }
    fprintf(gp, "e\n");
    fprintf(gp, "unset multiplot\n");
    fclose(gp);
    gp = popen("gnuplot -peristent", "w");
    if (!gp){
        cerr<<"Error: Could not open gp."<<endl;
        return 1;
    }
    fprintf(gp, "set terminal x11\n");
    fprintf(gp, "set grid\n");
    mu = 0.1;
    // animating mu goes from 0.1 to 4
    for (int frame=0; frame<100; frame++){
        fprintf(gp, "plot '-' with lines title 'RK4 method'\n");
        x_rk.clear();
        y_rk.clear();
        rk4Method(x_rk, y_rk, h, mu, t, x, y, steps);
        for (int i=0; i<steps; i++){
            fprintf(gp, "%f %f\n", x_rk[i], y_rk[i]);
        }
        fprintf(gp, "e\n");
        fflush(gp);
        usleep(100000);
        mu += 0.04;
    }
    fclose(gp);
    gp = popen("gnuplot -peristent", "w");
    if (!gp){
        cerr<<"Error: Could not open gp."<<endl;
        return 1;
    }
    fprintf(gp, "set terminal x11\n");
    fprintf(gp, "set grid\n");
    mu = 7;
    float x_sample[5] = {0, 2, 1.5, -1.5, 3}; //random values taken for initial conditions of x and y
    float y_sample[5] = {1, 0, -0.5, 0.5, 0};
    for (int frame=0; frame<5; frame++){
        fprintf(gp, "plot '-' with lines title 'Random initial conditions'\n");
        x = x_sample[frame];
        y = y_sample[frame];
        x_rk.clear();
        y_rk.clear();
        rk4Method(x_rk, y_rk, h, mu, t, x, y, steps);
        for (int i=0; i<steps; i++){
            fprintf(gp, "%f %f\n", x_rk[i], y_rk[i]);
        }
        fprintf(gp, "e\n");
        fflush(gp);
        usleep(2000000);
    }
    fclose(gp);
    
}




