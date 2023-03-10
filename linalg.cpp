#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>

using namespace std;

// TODO: LU, redo OOP (1.0), parallelism (2.0), arrays and pointers (3.0)

namespace linalg{

vector<vector<double>> rcdrop(vector<vector<double>> A, int r, int c){
    // get matrix by deleting a row and column, used for determinants and inverses
    vector<vector<double>> new_mat(c-1, vector<double>(r-1));
    int adj_i, adj_j;

    for(int i = 0; i < A.size(); i++){
        if(i != r){
            for(int j = 0; j < A[0].size(); j++){
                if(i != r && j != c){
                    if(i > r){adj_i = i-1;}else{adj_i = i;}
                    if(j > c){adj_j = j-1;}else{adj_j = j;}
                    new_mat[adj_i][adj_j] = A[i][j];
                }
            }
        } 
    }  
    return new_mat;
}

vector<vector<double>> matmul(vector<vector<double>> A, vector<vector<double>> B){
    // Matrix multimplication
    if(A[0].size() != B.size()){
        throw "num A rows must == num B cols";
    }
    vector<vector<double>> new_mat(A.size(), vector<double>(B[0].size()));
    for(int row = 0; row < A.size(); row++){
        for(int col = 0; col < B[0].size(); col++){
            double cur_val = 0;
            for(int ind = 0; ind < A[0].size(); ind ++){
                cur_val += A[row][ind] * B[ind][col]; // Matrix i,j value is the dot product of A row and B column
            }
            new_mat[row][col] = cur_val;
        }

    }
    return new_mat;
}

vector<vector<double>> scalar_mult(vector<vector<double>> A, double c){
    // multiply each element in matrix by a scalar
    for(int i = 0; i < A.size(); i++){
        for(int l = 0; l < A[0].size(); l++){
            A[i][l] = A[i][l] * c;
        }
    }
    return A;
}

vector<vector<double>> transpose(vector<vector<double>> A){

    //Swap rows and columns of matrix
    
    int rows = A[0].size();
    int cols = A.size();

    vector<vector<double>> new_mat(A[0].size(), vector<double>(A.size()));

    for(int r = 0; r<rows; r++){
        for(int c = 0; c<cols; c++){
            new_mat[r][c] = A[c][r];
        }
    }
    return new_mat;
}

void print_matrix(vector<vector<double>> A){
    // Output matrix MxN
    for(int i = 0; i < A.size(); i++){
        for(int l = 0; l < A[0].size(); l++){
            cout << A[i][l] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

double get_determinant(vector<vector<double>> A){
    if(A[0].size() != A.size()){
        throw "A must be square NxN";
    }
    if(A.size() == 1 && A[0].size() == 1){
        return A[0][0];
    }
    // Might not need this section, test
    if(A.size() == 2 && A[0].size() == 2){
        double det = (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
        return det;
    }
    double det = 0.0;
    int sign = 1;

    for(int l = 0; l < A[0].size(); l++){
        vector<vector<double>> nm = rcdrop(A, 0, l);
        det += sign*A[0][l]*get_determinant(nm); // alternate signs and multiply by determinant from dropped r,c pair
        sign = sign * -1;
    }
    return det;
}


vector<vector<double>> inverse(vector<vector<double>> A){
    
    if(A.size() == 2 && A[0].size() == 2){
        double det = (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
        if (det == 0){
            throw "Determinant is 0, no inverse exists";
        }
        return scalar_mult({{A[1][1], -1*A[0][1]}, {-1*A[1][0], A[0][0]}}, 1/det); // Inverse is 1/determinant scalar mult matrix
    }

    double det = get_determinant(A);
    if (det == 0){
            throw "Determinant is 0, no inverse exists";
    }
    vector<vector<double>> newMat = A;
    int sign = 1;
    for(int i=0; i<A.size(); i++){
        for(int l=0; l<A[0].size(); l++){
            if((i+1) % 2 == (l+1) % 2){
                sign = 1; // accomplish alternating signs without being affected by transpose below
            }
            else{
                sign = -1;
            }
            vector<vector<double>> nm = rcdrop(A, i, l);
            newMat[l][i] = get_determinant(nm) * sign * (1/det); // accomplish transpose by setting index [l][i]
        }
    }
    return newMat;
}

int min(int a, int b){
    if(a <= b){return a;}
    else{return b;}
}

vector<double> diagonal(vector<vector<double>> A){
    // switch to min of rown and cols
    int len = min(A.size(), A[0].size());
    vector<double> diag(len);
    for(int i=0; i < len; i++){
        diag[i] = A[i][i];
    }
    return diag;
}

int factorial(int x){
    int res = x;
    for(int i=2; i < x; ++i){
        res = res*i;
    }
    return res;
}

vector<double> trivec(vector<vector<double>> A){
    // gets all values of lower triangle
    vector<double> tri(factorial(min(A.size(), A[0].size()) - 1));
    for(int i = 0; i < A[0].size(); i++){
        for(int l = i+1; l < A.size(); l++){
            tri[i] = A[l][i];
        }
    }
    return tri;
}

bool trivec_zero(vector<double> tri, double tol){
    // checks if all values of lower triangle are near zero
    bool lt = true;
    for(int i=0; i < tri.size(); i++){
        if (abs(tri[i]) >= tol){
            return false;
        }
    }
    return lt;
}

double magnitude(vector<double> v){
    double s;
    for(int i=0; i < v.size(); i++){
        s += pow(v[i],2);
    }
    return sqrt(s);
}

vector<double> vector_scalar_mult(vector<double> v, double c){
    for(int i=0; i < v.size(); i++){
        v[i] = v[i] * c;
    }
    return v;
}

vector<double> unit_vector(vector<double> v){
    // divide vector by its magnitude
    double len = magnitude(v);
    if(len == 0){
        throw "Length is 0, magnitude is invalid";
    }
    v = vector_scalar_mult(v, 1/len);
    return v;
}

double vector_dot(vector<double> u, vector<double> v){
    // dot product
    // check vectors are same size
    if(v.size() != u.size()){
        throw "Vectors must be the same size!";
    }
    double dp;
    for(int i=0; i < v.size(); i++){
        dp += v[i] * u[i];
    }
    return dp;
}


vector<double> projection(vector<double> u, vector<double> v){
    // vector projection as defines by graham-schmidt process
    double scal = vector_dot(v,u) / vector_dot(u,u);
    return vector_scalar_mult(u, scal);
}

vector<double> vector_subtraction(vector<double> u, vector<double> v){
    // check vectors are same size
    if(v.size() != u.size()){
        throw "Vectors must be the same size!";
    }
    for(int i = 0; i < u.size(); i++){
        u[i] = u[i] - v[i];
    }
    return u;
}

struct AQR_dec {
    vector<vector<double>> A, Q, R;
};

AQR_dec QRDecomposition(vector<vector<double>> A){
    // QR Decomposition and return hash map of A, Q, and R
    // check if A is square
    if(A[0].size() != A.size()){
        throw "A must be square NxN";
    }

    vector<vector<double>> Q(A.size(), vector<double>(A[0].size())); // size dims of transpose
    vector<vector<double>> R, Qt;

    vector<vector<double>> At = transpose(A);
    vector<vector<double>> u(At[0].size(), vector<double> (At.size()));

    for(int i=0; i < At.size(); i++){
        u[i] = At[i]; // Use transpose to get A's column vectors easier
        for(int l=0; l < i; l++){
            u[i] = vector_subtraction(u[i], projection(u[l], At[i]));
        }
        Q[i] = unit_vector(u[i]);
    }

    AQR_dec AQR;

    Qt = Q; // since we were using the transpose to get column vectors
    Q = transpose(Q);
    u = transpose(u);
    R = matmul(Qt, A);
    A = matmul(R, Q);

    AQR.A = A;
    AQR.Q = Q;
    AQR.R = R;

    return AQR;
}

AQR_dec QRAlgorithm(vector<vector<double>> A, const int max_iter){
    // Run QR decomposition until A is upper-triangular
    int iter = 1;
    AQR_dec AQRi;
    vector<double> tv = trivec(A);
    bool tv0 = trivec_zero(tv, 1e-10); 

    while(iter <= max_iter && tv0 == false){
        AQRi = QRDecomposition(A);
        A = AQRi.A;
        tv = trivec(A);
        tv0 = trivec_zero(tv, 1e-10);
        iter += 1;
    }
    if(iter >= max_iter){
        cout << "\nWARNING: QRA did not converge! Set max_iter higher or use Wilkinson shift: A -> A - cI";
    }

    cout << "\n\nFINISHED QRA";
    return AQRi;
}

void print_vec(vector<double> v){
    for(int i=0; i < v.size(); i++){
        cout << "\n" << v[i];
    }
}

vector<double> eigenvalues(vector<vector<double>> A){
    // Eigenvalues are diagonal of A after running QR algorithm
    AQR_dec AQR = QRAlgorithm(A, 1000);
    A = AQR.A;
    return diagonal(A);
}

vector<vector<double>> cholesky(vector<vector<double>> A){
    // Can extend to check if positive definite
    if(A[0].size() != A.size()){
        throw "A must be square NxN";
    }
    vector<vector<double>> L = A;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++){
                sum += L[i][k] * L[j][k];
            }

            if (i == j){
                L[i][j] = sqrt(A[i][i] - sum);
            }
            else{
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
            }
        }
    }
    return L;
}

vector<vector<double>> diagonalize(vector<vector<double>> A){
    // P = matrix with eigenvalues as columns aka Q from QR algorithm
    // D = matrix with eigenvalues on diagonal
    vector<vector<double>> P, Pi, PDPi;
    vector<vector<double>> D(A[0].size(), vector<double>(A.size()));
    AQR_dec AQR = QRAlgorithm(A, 1000);
    P = AQR.Q;
    vector<double> eigenvals = diagonal(AQR.A);

    // Set the D matrix to diagonal matrix of eigenvals
    for(int i=0; i<AQR.A.size(); i++){
        for(int l=0; l<AQR.A.size(); l++){
            if(i != l){
                D[i][l] = 0;
            }
            else{
                D[i][l] = eigenvals[i];
            }
        }
    }

    Pi = inverse(P);
    PDPi = matmul(P,D);
    PDPi = matmul(PDPi, Pi);
    return D;
}

struct SVD_dec{
    vector<vector<double>> U, VT, Sigma;
};

SVD_dec SVDecomp(vector<vector<double>> A){
    // Run QR decomposition until A is upper-triangular
    vector<vector<double>> AT = linalg::transpose(A);
    vector<vector<double>> ATA = linalg::matmul(AT, A);
    vector<vector<double>> AAT = linalg::matmul(A, AT);
    
    AQR_dec ATA_AQR = QRAlgorithm(ATA, 1000);
    AQR_dec AAT_AQR = QRAlgorithm(AAT, 1000);
    vector<double> sigma2s = diagonal(ATA_AQR.A);
    vector<vector<double>> Sigma = A;

    for(int i = 0; i < A.size(); i++){
        for(int l = 0; l < A.size(); l++){
            if(i == l){
                Sigma[i][l] = sqrt(sigma2s[i]);
            }
            else{
                Sigma[i][l] = 0;
            }
        }
    }

    SVD_dec SVD;
    SVD.U = AAT_AQR.Q;
    SVD.VT = transpose(ATA_AQR.Q);
    SVD.Sigma = Sigma;
    return SVD;
}
}


class Matrix{
    // Manipulate as a class 
    public:
    vector<vector<double>> mat;
    vector<double> eigenvals;
    vector<vector<double>> eigenvectors;

    void set_mat(vector<vector<double>> A){
        mat = A;
    }

    void matmul(Matrix B){
        mat = linalg::matmul(mat, B.mat);
    }

    void scalar_mult(double c){
        mat = linalg::scalar_mult(mat, c);
    }

    void transpose(){
        mat = linalg::transpose(mat);
    }

    void inverse(){
        mat = linalg::inverse(mat);
    }

    void get_eigens(){
        linalg::AQR_dec AQR = linalg::QRAlgorithm(mat, 1000);
        eigenvals = linalg::diagonal(AQR.A);
        eigenvectors = AQR.Q;
    }

    void print_matrix(){
        linalg::print_matrix(mat);
    }

};


int main(){
    // MATMUL TESTS

    linalg::print_matrix(linalg::matmul({{1}}, {{7}}));
    linalg::print_matrix(linalg::matmul({{1,2,3}, {4,5,6}}, {{10, 11}, {20,21}, {30,31}}));
    linalg::print_matrix(linalg::matmul({{1,2,3,4}}, {{10, 11}, {20,21}, {30,31}, {4,5}}));

    // DETERMINANT TESTS
    cout << "\nDeterminant: " << linalg::get_determinant({{-5, -4}, {-2,-3}});
    cout << "\nDeterminant: " << linalg::get_determinant({{2,-3,1}, {2,0,-1}, {1,4,5}});
    cout << "\nDeterminant: " << linalg::get_determinant({{4,3,2,2}, {0,1,-3,3}, {0,-1,3,3}, {0,3,1,1}});

    // INVERSE TESTS
    cout << "\n\n";
    linalg::print_matrix(linalg::inverse({{4,7}, {2,6}}));
    linalg::print_matrix(linalg::inverse({{2,-3,1}, {2,0,-1}, {1,4,5}}));
    linalg::print_matrix(linalg::inverse({{2,6,4,-1}, {3,0,7,4},{-3,3,4,1}, {8,-1,1,0}}));
    linalg::print_matrix(linalg::inverse({{-3,-1,2,-3}, {-3,1,2,-2},{-2,3,0,1}, {1,-2,-3,1}}));
    linalg::print_matrix(linalg::inverse({{2,3,7,-1,-1}, {1,4,-2,0,1},{3,2,2,-1,3}, {0,1,0,4,7}, {5,2,-3,-1,1}}));

    // EIGENVALUE TESTS
    cout << "\nEIGENVALS";
    linalg::print_vec(linalg::eigenvalues({{1,2,3}, {3,2,1}, {2,1,3}}));
    linalg::print_vec(linalg::eigenvalues({{4,8,1},{0,1,0},{2,-3,-1}}));
    linalg::print_vec(linalg::eigenvalues({{3,1,2,4}, {0,1,0,-2}, {5,2,2,2}, {3,4,0,1}}));

    linalg::AQR_dec AQR = linalg::QRAlgorithm({{1,2,3}, {3,2,1}, {2,1,3}}, 10000);
    linalg::print_matrix(AQR.A);

    // SVD
    
    cout << "\nSVD";
    linalg::SVD_dec SVD = linalg::SVDecomp({{1,2,3}, {3,2,1}, {2,1,3}, {-1,2,0}});
    linalg::print_matrix(SVD.Sigma);

    return 0;
}
