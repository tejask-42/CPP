#include <iostream>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <unistd.h> 
#include <tuple>
using namespace std;

float rho=28, sigma=10, b=float(8)/3;
float x_dot(float sigma, float x, float y){ return sigma*(y - x);}
float y_dot(float rho, float x, float y, float z){ return x*(rho-z) - y;}
float z_dot(float b, float x, float y, float z){ return x*y - b*z;}

// Euler Method
void eulerMethod(vector<float>& x_eul, vector<float>& y_eul, vector<float>& z_eul, float h, float t, float x, float y, float z, int steps, float rho, float sigma, float b) {
    for (int i = 0; i < steps; i++) {
        x_eul.push_back(x);
        y_eul.push_back(y);
        z_eul.push_back(z);

        x += h * x_dot(sigma, x, y);
        y += h * y_dot(rho, x, y, z);
        z += h * z_dot(b, x, y, z);
        t += h;
    }
}
// RK4 Method
void rk4Method(vector<float>& x_rk, vector<float>& y_rk, vector<float>& z_rk, float h, float t, float x, float y, float z, int steps, float rho, float sigma, float b) {
    for (int i = 0; i < steps; i++) {
        x_rk.push_back(x);
        y_rk.push_back(y);
        z_rk.push_back(z);

        float k1x = x_dot(sigma, x, y);
        float k1y = y_dot(rho, x, y, z);
        float k1z = z_dot(b, x, y, z);

        float x_k2 = x + h * k1x / 2;
        float y_k2 = y + h * k1y / 2;
        float z_k2 = z + h * k1z / 2;
        float k2x = x_dot(sigma, x_k2, y_k2);
        float k2y = y_dot(rho, x_k2, y_k2, z_k2);
        float k2z = z_dot(b, x_k2, y_k2, z_k2);

        float x_k3 = x + h * k2x / 2;
        float y_k3 = y + h * k2y / 2;
        float z_k3 = z + h * k2z / 2;
        float k3x = x_dot(sigma, x_k3, y_k3);
        float k3y = y_dot(rho, x_k3, y_k3, z_k3);
        float k3z = z_dot(b, x_k3, y_k3, z_k3);

        float x_k4 = x + h * k3x;
        float y_k4 = y + h * k3y;
        float z_k4 = z + h * k3z;
        float k4x = x_dot(sigma, x_k4, y_k4);
        float k4y = y_dot(rho, x_k4, y_k4, z_k4);
        float k4z = z_dot(b, x_k4, y_k4, z_k4);

        x += h * (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y += h * (k1y + 2 * k2y + 2 * k3y + k4y) / 6;

        t += h;
    }
}

int main(){
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        cerr << "Error opening pipe to gp" << endl;
        return 1;
    }

    float h = 0.01, t=0, x=1, y=1, z=1;
    int steps = 1000;
    vector<float> x_eul, y_eul, z_eul, x_rk, y_rk, z_rk, time;
    eulerMethod(x_eul, y_eul, z_eul, h, t, x, y, z, steps, rho, sigma, b);
    rk4Method(x_rk, y_rk, z_rk, h, t, x, y, z, steps, rho, sigma, b);

    float t_ind = 0;
    for (int i=0; i<steps; i++){
        time.push_back(t_ind);
        t_ind += h;
    }
    fprintf(gp, "set multiplot layout 2,3 rowsfirst\n");

    fprintf(gp, "set title 'Lorenz Attractor: Euler vs RK4'\n");
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

    fprintf(gp, "set title 'Lorenz Attractor: Euler vs RK4'\n");
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

    fprintf(gp, "set title 'Lorenz Attractor: Euler vs RK4'\n");
    fprintf(gp, "set xlabel 'Time'\n");
    fprintf(gp, "set ylabel 'z'\n");
    fprintf(gp, "plot '-' with lines title 'Euler z', '-' with lines title 'RK4 z'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], z_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", time[i], z_rk[i]);
    }
    fprintf(gp, "e\n");
    
    fprintf(gp, "set title 'y vs x for rho=28'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "plot '-' with lines title 'Euler', '-' with lines title 'RK4'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", x_eul[i], y_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < 2 * steps; i++) {
        fprintf(gp, "%f %f\n", x_rk[i], y_rk[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set title 'z vs y for rho=28'\n");
    fprintf(gp, "set xlabel 'y'\n");
    fprintf(gp, "set ylabel 'z'\n");
    fprintf(gp, "plot '-' with lines title 'Euler', '-' with lines title 'RK4'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", y_eul[i], z_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < 2 * steps; i++) {
        fprintf(gp, "%f %f\n", y_rk[i], z_rk[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set title 'x vs z for rho=28'\n");
    fprintf(gp, "set xlabel 'z'\n");
    fprintf(gp, "set ylabel 'x'\n");
    fprintf(gp, "plot '-' with lines title 'Euler', '-' with lines title 'RK4'\n");
    for (int i = 0; i < steps; i++) {
        fprintf(gp, "%f %f\n", z_eul[i], x_eul[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < 2 * steps; i++) {
        fprintf(gp, "%f %f\n", z_rk[i], x_rk[i]);
    }
    fprintf(gp, "e\n");
    fprintf(gp, "unset multiplot\n");
    pclose(gp);    
}
