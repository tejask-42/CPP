#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <eigen/Eigen/Sparse>
using namespace std;

void make_dist(vector<vector<double>> &grid, int n, int dist = 1){
    if (dist == 1){
        grid[n / 2 + 10][n / 2] = 1.0;
        grid[n / 2 - 10][n / 2] = -1.0;

        for (int j=0; j<n; j++){
            grid[0][j] = 1.0;
            grid[n-1][j] = -1.0;
        }
    }
    else if (dist == 2){ // quadrupole
        grid[n / 2 + 10][n / 2 + 10] = 1.0;
        grid[n / 2 + 10][n / 2 - 10] = -1.0;
        grid[n / 2 - 10][n / 2 + 10] = -1.0;
        grid[n / 2 - 10][n / 2 - 10] = 1.0;
    }
    else if (dist == 3){ // ring of charge
        int radius = n / 4, cx = n / 2, cy = n / 2;
        int num = 200;
        for (int k=0; k < num; k++){
            double angle = 2 * M_PI * k / num;
            int x = cx + int(radius * cos(angle));
            int y = cy + int(radius * sin(angle));
            if (x >= 0 && x < n && y >= 0 && y < n){
                grid[x][y] = 1.0;
            }
        }
    }
}

void flatten(const vector<vector<double>> &grid, Eigen::VectorXd &flat){
    int n = grid.size();
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            flat[i * n + j] = grid[i][j];
        }
    }
}

void flat2grid(const Eigen::VectorXd &flat, vector<vector<double>> &grid){
    int n = sqrt(flat.size());
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            grid[i][j] = flat[i * n + j];
        }
    }
}

void make_laplacian_operator(Eigen::SparseMatrix<double>& op_matrix, int n) {
    int size = n * n; 
    std::vector<Eigen::Triplet<double>> triplets; 

    double delta = -1.0 / ((n - 1) * (n - 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = i * n + j;
            triplets.emplace_back(row, row, -4.0 * delta);
            if (i + 1 < n) triplets.emplace_back(row, (i + 1) * n + j, delta); 
            if (i - 1 >= 0) triplets.emplace_back(row, (i - 1) * n + j, delta); 
            if (j + 1 < n) triplets.emplace_back(row, i * n + (j + 1), delta); 
            if (j - 1 >= 0) triplets.emplace_back(row, i * n + (j - 1), delta); 
        }
    }

    // Build the sparse matrix from triplets
    op_matrix.setFromTriplets(triplets.begin(), triplets.end());
}

void solve_system_lineq(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, int maxIter = 200, double tol = 1e-5) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.setMaxIterations(maxIter);
    cg.setTolerance(tol);
    cg.compute(A);
    x = cg.solve(b);
    
    if (cg.info() == Eigen::Success) {
        std::cout << "Converged after " << cg.iterations() << " iterations." << std::endl;
    } else {
        std::cout << "Conjugate gradient did not converge." << std::endl;
    }
}

int main(){
    int n = 50;
    vector<vector<double>> grid(n, vector<double>(n, 0.0));
    Eigen::VectorXd p_flat(n * n);
    make_dist(grid, n, 3); // defining charge distribution
    flatten(grid, p_flat); // flattening it to get p_flat
    Eigen::SparseMatrix<double> op_matrix(n * n, n * n);
    make_laplacian_operator(op_matrix, n); // M matrix
    Eigen::VectorXd u_flat = Eigen::VectorXd::Zero(n * n); 
    solve_system_lineq(op_matrix, p_flat, u_flat, n * n); 
    vector<vector<double>> u_grid(n, vector<double>(n, 0.0));
    flat2grid(u_flat, u_grid);
    std::ofstream field_file("field.txt");
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            field_file << u_grid[i][j] << " "; 
        }
        field_file << "\n";
    }
    field_file.close();

}
