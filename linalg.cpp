#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>

using namespace std;

namespace linalg
{

    vector<vector<double>> rcdrop(vector<vector<double>> &A, int r, int c)
    {
        // get matrix by deleting a row and column, used for determinants and inverses
        vector<vector<double>> newMat(A.size() - 1, vector<double>(A[0].size() - 1, 0)); // why on earth does not initializing 0 cause a bug?
        int adj_i, adj_j;
        for (int i = 0; i < A.size(); i++)
        {
            if (i != r)
            {
                for (int j = 0; j < A[0].size(); j++)
                {
                    if (j != c)
                    {
                        if (i > r)
                        {
                            adj_i = i - 1;
                        }
                        else
                        {
                            adj_i = i;
                        }
                        if (j > c)
                        {
                            adj_j = j - 1;
                        }
                        else
                        {
                            adj_j = j;
                        }
                        newMat[adj_i][adj_j] = A[i][j];
                    }
                }
            }
        }
        return newMat;
    }

    vector<vector<double>> matmul(vector<vector<double>> &A, vector<vector<double>> &B)
    {
        // Matrix multimplication
        if (A[0].size() != B.size())
        {
            throw "num A rows must == num B cols";
        }
        vector<vector<double>> newMat(A.size(), vector<double>(B[0].size(), 0));
        for (int row = 0; row < A.size(); row++)
        {
            for (int col = 0; col < B[0].size(); col++)
            {
                double curVal = 0;
                for (int ind = 0; ind < A[0].size(); ind++)
                {
                    curVal += A[row][ind] * B[ind][col]; // Matrix i,j value is the dot product of A row and B column
                }
                newMat[row][col] = curVal;
            }
        }
        return newMat;
    }

    vector<vector<double>> scalar_mult(vector<vector<double>> A, double c)
    {
        // multiply each element in matrix by a scalar
        for (int i = 0; i < A.size(); i++)
        {
            for (int l = 0; l < A[0].size(); l++)
            {
                A[i][l] = A[i][l] * c;
            }
        }
        return A;
    }

    vector<vector<double>> transpose(vector<vector<double>> &A)
    {

        // Swap rows and columns of matrix

        int rows = A[0].size();
        int cols = A.size();

        vector<vector<double>> newMat(rows, vector<double>(cols)); // AAAGGGHH Why does initializing with 0 cause a bug here?

        for (int r = 0; r < rows; r++)
        {
            for (int c = 0; c < cols; c++)
            {
                newMat[r][c] = A[c][r];
            }
        }
        return newMat;
    }

    void print_matrix(vector<vector<double>> &A)
    {
        // Output matrix MxN
        for (int i = 0; i < A.size(); i++)
        {
            for (int l = 0; l < A[0].size(); l++)
            {
                cout << A[i][l] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }

    double get_determinant(vector<vector<double>> &A)
    {
        if (A[0].size() != A.size())
        {
            throw "A must be square NxN";
        }
        if (A.size() == 1 && A[0].size() == 1)
        {
            return A[0][0];
        }
        // Might not need this section, test
        if (A.size() == 2 && A[0].size() == 2)
        {
            double det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
            return det;
        }
        double det = 0.0;
        int sign = 1;

        for (int l = 0; l < A[0].size(); l++)
        {
            vector<vector<double>> nm = rcdrop(A, 0, l);
            det += sign * A[0][l] * get_determinant(nm); // alternate signs and multiply by determinant from dropped r,c pair
            sign = sign * -1;
        }

        return det;
    }

    vector<vector<double>> inverse(vector<vector<double>> &A)
    {

        if (A.size() == 2 && A[0].size() == 2)
        {
            double det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
            if (det == 0)
            {
                throw "Determinant is 0, no inverse exists";
            }
            vector<vector<double>> Adet = {{A[1][1], -1 * A[0][1]}, {-1 * A[1][0], A[0][0]}};
            return scalar_mult(Adet, 1 / det); // Inverse is 1/determinant scalar mult matrix
        }

        double det = get_determinant(A);
        if (det == 0)
        {
            throw "Determinant is 0, no inverse exists";
        }

        vector<vector<double>> newMat = A;
        int sign = 1;
        for (int i = 0; i < A.size(); i++)
        {
            for (int l = 0; l < A[0].size(); l++)
            {
                if ((i + 1) % 2 == (l + 1) % 2)
                {
                    sign = 1; // accomplish alternating signs without being affected by transpose below
                }
                else
                {
                    sign = -1;
                }
                vector<vector<double>> nm = rcdrop(A, i, l);
                newMat[l][i] = get_determinant(nm) * sign * (1 / det); // accomplish transpose by setting index [l][i]
            }
        }
        return newMat;
    }

    int min(int a, int b)
    {
        if (a >= b)
        {
            return a;
        }
        else
        {
            return b;
        }
    }

    vector<double> diagonal(vector<vector<double>> &A)
    {
        vector<double> diag(min(A[0].size(), A.size()));
        for (int i = 0; i < A.size(); i++)
        {
            diag[i] = A[i][i];
        }
        return diag;
    }

    int sum_digits(int n)
    { // inclusive
        return (n * (n + 1)) / 2;
    }

    vector<double> trivec(vector<vector<double>> &A)
    {
        // gets all values of lower triangle
        vector<double> tri(sum_digits(min(A.size(), A[0].size()) - 1));
        int counter = 0;
        for (int i = 0; i < A.size(); i++)
        {
            for (int l = 0; l < i; l++)
            {
                tri[counter] = A[i][l];
                counter++;
            }
        }
        return tri;
    }

    bool trivec_zero(vector<double> &tri, double tol)
    {
        // checks if all values of lower triangle are near zero
        bool lt = true;
        for (int i = 0; i < tri.size(); i++)
        {
            if (abs(tri[i]) >= tol)
            {
                return false;
            }
        }
        return lt;
    }

    double magnitude(vector<double> &v)
    {
        double s = 0;
        for (int i = 0; i < v.size(); i++)
        {
            s += pow(v[i], 2);
        }
        return sqrt(s);
    }

    vector<double> vector_scalar_mult(vector<double> v, double c)
    {
        for (int i = 0; i < v.size(); i++)
        {
            v[i] = v[i] * c;
        }
        return v;
    }

    vector<double> unit_vector(vector<double> v)
    {
        // divide vector by its magnitude
        double len = magnitude(v);
        if (len == 0)
        {
            throw "Length is 0, magnitude is invalid";
        }
        v = vector_scalar_mult(v, 1 / len);
        return v;
    }

    double vector_dot(vector<double> u, vector<double> v)
    {
        // dot product
        // check vectors are same size
        if (v.size() != u.size())
        {
            throw "Vectors must be the same size!";
        }
        double dp = 0;
        for (int i = 0; i < v.size(); i++)
        {
            dp += v[i] * u[i];
        }
        return dp;
    }

    vector<double> projection(vector<double> u, vector<double> v)
    {
        // vector projection as defines by graham-schmidt process
        double scal = vector_dot(v, u) / vector_dot(u, u);
        return vector_scalar_mult(u, scal);
    }

    vector<double> vector_subtraction(vector<double> u, vector<double> v)
    {
        // check vectors are same size
        if (v.size() != u.size())
        {
            throw "Vectors must be the same size!";
        }
        for (int i = 0; i < u.size(); i++)
        {
            u[i] = u[i] - v[i];
        }
        return u;
    }

    unordered_map<string, vector<vector<double>>> QRDecomposition(vector<vector<double>> A)
    {
        // QR Decomposition and return hash map of A, Q, and R
        // check if A is square
        if (A[0].size() != A.size())
        {
            throw "A must be square NxN";
        }

        vector<vector<double>> Q(A.size(), vector<double>(A[0].size(), 0));
        vector<vector<double>> u(A.size(), vector<double>(A[0].size(), 0));
        vector<vector<double>> R, Qt;
        vector<vector<double>> At = transpose(A);

        for (int i = 0; i < At.size(); i++)
        {
            vector<double> ui = At[i]; // Use transpose to get A's column vectors easier
            for (int l = 0; l < u.size(); l++)
            {
                ui = vector_subtraction(ui, projection(u[l], At[i]));
            }
            u[i] = ui;
            Q[i] = unit_vector(ui);
        }
        unordered_map<string, vector<vector<double>>> AQR;

        Qt = Q; // since we were using the transpose to get column vectors
        Q = transpose(Q);
        u = transpose(u);
        R = matmul(Qt, A);
        A = matmul(R, Q);

        AQR["A"] = A;
        AQR["Q"] = Q;
        AQR["R"] = R;

        return AQR;
    }

    unordered_map<string, vector<vector<double>>> QRAlgorithm(vector<vector<double>> A, const int max_iter)
    {
        // Run QR decomposition until A is upper-triangular
        int iter = 1;
        unordered_map<string, vector<vector<double>>> AQRi;
        vector<double> tv = trivec(A);
        bool tv0 = trivec_zero(tv, 1e-10);

        while (iter <= max_iter && tv0 == false)
        {
            AQRi = QRDecomposition(A);
            A = AQRi["A"];
            tv = trivec(A);
            tv0 = trivec_zero(tv, 1e-10);
            iter += 1;
        }
        if (iter >= max_iter)
        {
            cout << "\nWARNING: QRA did not converge! Set max_iter higher or use Wilkinson shift: A -> A - cI";
        }

        cout << "\n\nFINISHED QRA";
        return AQRi;
    }

    void print_vec(vector<double> &v)
    {
        for (int i = 0; i < v.size(); i++)
        {
            cout << "\n"
                 << v[i];
        }
    }

    vector<double> eigenvalues(vector<vector<double>> A)
    {
        // Eigenvalues are diagonal of A after running QR algorithm
        unordered_map<string, vector<vector<double>>> AQR = QRAlgorithm(A, 1000);
        return diagonal(AQR["A"]);
    }

    vector<vector<double>> cholesky(vector<vector<double>> A)
    {
        // Can extend to check if positive definite
        if (A[0].size() != A.size())
        {
            throw "A must be square NxN";
        }
        vector<vector<double>> L = A;
        for (int i = 0; i < A.size(); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += L[i][k] * L[j][k];
                }

                if (i == j)
                {
                    L[i][j] = sqrt(A[i][i] - sum);
                }
                else
                {
                    L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
                }
            }
        }
        return L;
    }

    vector<vector<double>> diagonalize(vector<vector<double>> A)
    {
        // P = matrix with eigenvectors as columns aka Q from QR algorithm
        // D = matrix with eigenvalues on diagonal
        vector<vector<double>> P, Pi, PDPi;
        unordered_map<string, vector<vector<double>>> AQR = QRAlgorithm(A, 1000);
        vector<vector<double>> D(AQR["A"].size(), vector<double>(AQR["A"].size()));
        P = AQR["Q"];
        vector<double> eigenvals = diagonal(AQR["A"]);

        // Set the D matrix to diagonal matrix of eigenvals
        for (int i = 0; i < AQR["A"].size(); i++)
        {
            for (int l = 0; l < AQR["A"].size(); l++)
            {
                if (i != l)
                {
                    D[i][l] = 0;
                }
                else
                {
                    D[i][l] = eigenvals[i];
                }
            }
        }

        Pi = inverse(P);
        PDPi = matmul(P, D);
        PDPi = matmul(PDPi, Pi);
        return D;
    }

    unordered_map<string, vector<vector<double>>> SVDecomp(vector<vector<double>> A)
    {
        // Run QR decomposition until A is upper-triangular
        vector<vector<double>> AT = linalg::transpose(A);
        vector<vector<double>> ATA = linalg::matmul(AT, A);
        vector<vector<double>> AAT = linalg::matmul(A, AT);

        unordered_map<string, vector<vector<double>>> ATA_AQR = QRAlgorithm(ATA, 1000);
        unordered_map<string, vector<vector<double>>> AAT_AQR = QRAlgorithm(AAT, 1000);
        vector<double> sigma2s = diagonal(ATA_AQR["A"]);
        vector<vector<double>> Sigma = A;

        for (int i = 0; i < A.size(); i++)
        {
            for (int l = 0; l < A.size(); l++)
            {
                if (i == l)
                {
                    Sigma[i][l] = sqrt(sigma2s[i]);
                }
                else
                {
                    Sigma[i][l] = 0;
                }
            }
        }

        unordered_map<string, vector<vector<double>>> SVD;
        SVD["U"] = AAT_AQR["Q"];
        SVD["VT"] = transpose(ATA_AQR["Q"]);
        SVD["Sigma"] = Sigma;

        return SVD;
    }

}

class Matrix
{
    // Manipulate as a class
public:
    vector<vector<double>> mat;
    vector<double> eigenvals;
    vector<vector<double>> eigenvectors;

    void set_mat(vector<vector<double>> A)
    {
        mat = A;
    }

    void matmul(Matrix B)
    {
        mat = linalg::matmul(mat, B.mat);
    }

    void scalar_mult(double c)
    {
        mat = linalg::scalar_mult(mat, c);
    }

    void transpose()
    {
        mat = linalg::transpose(mat);
    }

    void inverse()
    {
        mat = linalg::inverse(mat);
    }

    void get_eigens()
    {
        unordered_map<string, vector<vector<double>>> AQR = linalg::QRAlgorithm(mat, 1000);
        eigenvals = linalg::diagonal(AQR["A"]);
        eigenvectors = AQR["Q"];
    }

    void print_matrix()
    {
        linalg::print_matrix(mat);
    }
};

int main()
{
    // MATMUL TESTS
    vector<vector<double>> a1 = {{1}};
    vector<vector<double>> a2 = {{7}};
    vector<vector<double>> am1 = linalg::matmul(a1, a2);
    linalg::print_matrix(am1);

    vector<vector<double>> a3 = {{1, 2, 3}, {4, 5, 6}};
    vector<vector<double>> a4 = {{10, 11}, {20, 21}, {30, 31}};
    vector<vector<double>> am2 = linalg::matmul(a3, a4);
    linalg::print_matrix(am2);

    vector<vector<double>> a5 = {{1, 2, 3, 4}};
    vector<vector<double>> a6 = {{10, 11}, {20, 21}, {30, 31}, {4, 5}};
    vector<vector<double>> am3 = linalg::matmul(a5, a6);
    linalg::print_matrix(am3);

    // DETERMINANT TESTS
    vector<vector<double>> a7 = {{-5, -4}, {-2, -3}};
    cout << "\nDeterminant: " << linalg::get_determinant(a7);

    vector<vector<double>> a8 = {{2, -3, 1}, {2, 0, -1}, {1, 4, 5}};
    cout << "\nDeterminant: " << linalg::get_determinant(a8);

    vector<vector<double>> a9 = {{4, 3, 2, 2}, {0, 1, -3, 3}, {0, -1, 3, 3}, {0, 3, 1, 1}};
    cout << "\nDeterminant: " << linalg::get_determinant(a9);

    // INVERSE TESTS
    cout << "\n\n";
    vector<vector<double>> a10 = {{4, 7}, {2, 6}};
    vector<vector<double>> i10 = linalg::inverse(a10);
    linalg::print_matrix(i10);

    vector<vector<double>> a11 = {{2, -3, 1}, {2, 0, -1}, {1, 4, 5}};
    vector<vector<double>> i11 = linalg::inverse(a11);
    linalg::print_matrix(i11);

    vector<vector<double>> a12 = {{2, 6, 4, -1}, {3, 0, 7, 4}, {-3, 3, 4, 1}, {8, -1, 1, 0}};
    vector<vector<double>> i12 = linalg::inverse(a12);
    linalg::print_matrix(i12);

    vector<vector<double>> a13 = {{-3, -1, 2, -3}, {-3, 1, 2, -2}, {-2, 3, 0, 1}, {1, -2, -3, 1}};
    vector<vector<double>> i13 = linalg::inverse(a13);
    linalg::print_matrix(i13);

    vector<vector<double>> a14 = {{2, 3, 7, -1, -1}, {1, 4, -2, 0, 1}, {3, 2, 2, -1, 3}, {0, 1, 0, 4, 7}, {5, 2, -3, -1, 1}};
    vector<vector<double>> i14 = linalg::inverse(a14);
    linalg::print_matrix(i14);

    // EIGENVALUE TESTS
    cout << "\nEIGENVALS";

    vector<vector<double>> a15 = {{1, 2, 3}, {3, 2, 1}, {2, 1, 3}};
    vector<double> e15 = linalg::eigenvalues(a15);
    linalg::print_vec(e15);

    vector<vector<double>> a16 = {{4, 8, 1}, {0, 1, 0}, {2, -3, -1}};
    vector<double> e16 = linalg::eigenvalues(a16);
    linalg::print_vec(e16);

    vector<vector<double>> a17 = {{3, 1, 2, 4}, {0, 1, 0, -2}, {5, 2, 2, 2}, {3, 4, 0, 1}};
    vector<double> e17 = linalg::eigenvalues(a17);
    linalg::print_vec(e17);

    vector<vector<double>> a18 = {{1, 2, 3}, {3, 2, 1}, {2, 1, 3}};
    unordered_map<string, vector<vector<double>>> AQR = linalg::QRAlgorithm(a18, 10000);
    linalg::print_matrix(AQR["A"]);

    // SVD

    cout << "\nSVD";
    unordered_map<string, vector<vector<double>>> SVD = linalg::SVDecomp({{1, 2, 3}, {3, 2, 1}, {2, 1, 3}, {-1, 2, 0}});
    vector<vector<double>> sig = SVD["Sigma"];
    linalg::print_matrix(sig);

    return 0;
}
