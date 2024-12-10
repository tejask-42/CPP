#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdio>
#include<vector>
using namespace std;
/*
    |4 -1 -1|       |x|        |3 |
A = |-2 6  1|   x = |y|    b = |9 |
    |1  1  7|       |z|        |-6|
*/ 
// Jacobi method
void jacobi(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, int n, int maxIter = 100, double tol = 1e-5) {
    vector<double> newX(n, 0.0); // stores new values of soln in each iteration

    for (int iter = 1; iter <= maxIter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) { sum -= A[i][j] * x[j]; }
            }
            newX[i] = sum / A[i][i];
        }
        // convergence check
        double maxDiff = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiff = max(maxDiff, fabs(newX[i] - x[i]));
        }
        x = newX;
        if (maxDiff < tol) {
            cout << "Converged after " << iter << " iterations.\n";
            return;
        }
    }
    cout << "Did not converge within the maximum number of iterations.\n";
}

// Gauss - Seidel method 
void gaussSeidel(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, int n, int maxIter = 100, double tol = 1e-5) {
    for (int iter = 1; iter <= maxIter; iter++) {
        vector<double> oldX = x; // storing previous iteration values

        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {sum -= A[i][j] * x[j];}
            }
            x[i] = sum / A[i][i]; 
        }
        // convergence check
        double maxDiff = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiff = max(maxDiff, fabs(x[i] - oldX[i]));
        }
        if (maxDiff < tol) {
            cout << "Converged after " << iter << " iterations.\n";
            return;
        }
    }
    cout << "Did not converge within the maximum number of iterations.\n";
}
int main() {
    vector<vector<double>> A(3, vector<double>(3));
    A = {
        {4, -1, -1},
        {-2, 6, 1},
        {-1, 1, 7}
    };

    vector<double> b = {3, 9, -6};

    // initial guess (0, 0, 0)
    vector<double> x(3, 0.0);
    cout<<"Applying Jacobi method:"<<endl;
    jacobi(A, b, x, 3);
    cout<<"Solution"<<endl;
    cout<<"x: "<<x[0]<<endl<<"y: "<<x[1]<<endl<<"z: "<<x[2]<<endl;
    cout<<"-----------------------------------------------"<<endl;
    x = {0.0, 0.0, 0.0}; // reinitialising initial guess to (0, 0, 0)
    cout<<"Applying Gauss-Seidel method:"<<endl;
    gaussSeidel(A, b, x, 3);

    cout << "Solution:"<<endl;
    cout<<"x: "<<x[0]<<endl<<"y: "<<x[1]<<endl<<"z: "<<x[2]<<endl;
   

}