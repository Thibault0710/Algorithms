#include <bits/stdc++.h>
using namespace std;

template <typename T>
void display(vector<vector<T>> &A) {
        for(int i = 0; i < A.size(); ++i) {
                for(int j = 0; j < A[0].size(); ++j) cout << A[i][j] << ' ';
                cout << endl;
        }
        cout << endl;
}

template <typename T>
void display(vector<T> &vec) {
        for(int i = 0; i < vec.size(); ++i) cout << vec[i] << ' ';
        cout << endl<<endl;
}

template <typename T>
void pivot(vector<vector<T>> &A, vector<T> &b, vector<T> &c, T &nu,int e, int l) {
        A[l][e] = 1/A[l][e];
        b[l]    = -b[l] * A[l][e];
        for(int i = 0; i < A[0].size(); ++i) if(i != e) A[l][i] *= -A[l][e];

        for(int i = 0; i < A.size(); ++i) {
                int line = i;
                if(line == l) continue;

                for(int j = 0; j < A[0].size(); ++j) if(j != e) A[i][j] += A[i][e]*A[l][j];
                b[i]    += A[i][e] * b[l];
                A[i][e]  *= A[l][e];

        }

        for(int i = 0; i < A[0].size(); ++i) if(i != e) c[i] += A[l][i] * c[e];
        nu   += c[e]*b[l];
        c[e] *= A[l][e];
}

template <typename T>
int continue_simplex(const vector<T> &c) {
        for(int i = 0; i < c.size(); ++i) {
                if(c[i] > 0) return i;
        }
        return -1;
}

template <typename T>
vector<T> simplex(vector<vector<T>> A, vector<T> b, vector<T> c) {
        for(int i = 0; i < A.size(); ++i) for(int j = 0; j < A[0].size(); ++j) A[i][j] = -A[i][j];

        vector<int> hbase(A[0].size());
        vector<int> base(A.size());
        for(int i = 0; i < A[0].size(); ++i)    hbase[i]  = i;
        for(int i = 0; i < A.size(); ++i) base[i] = A[0].size()+i;

        T nu = 0;
        int index = continue_simplex(c);
        while(index >= 0) {
                int eq = 0;
                T eq_value = (A[0][index] < 0) ?  -(b[0]/A[0][index]) : numeric_limits<T>::max();
                for(int i = 1; i < A.size(); ++i) {
                        if(A[i][index] >= 0) continue;
                        auto value = -(b[i] / A[i][index]);
                        if(value < eq_value) {
                                eq_value = value;
                                eq       = i;
                        }
                }
                if(eq_value == numeric_limits<T>::max()) {
                        cout <<"NON BORNE"<<endl;
                        break;
                }
                pivot(A,b,c,nu,index,eq);//index is the entering variable and eqs the outgoing one
                swap(hbase[index], base[eq]);
                index = continue_simplex(c);
        }

        vector<T> optimal(A[0].size());
        for(int i = 0; i < A.size(); ++i) if(base[i] < A[0].size()) optimal[base[i]] = b[i];

        cout << "Optimal solution: " << endl;
        display(optimal);
        cout << "Optimal value (nu) = " << nu << endl; // Print optimal value

        return optimal;
}

int main() {
        int m, n;
        cout << "Entrez le nombre de contraintes (m) : ";
        cin >> m;
        cout << "Entrez le nombre de variables (n) : ";
        cin >> n;

        vector<vector<double>> A(m, vector<double>(n));
        vector<double> b(m);
        vector<double> c(n);
        cout << "Entrez les valeurs de la matrice A (m x n) :" << endl;
        for(int i = 0; i < m; ++i) for(int j = 0; j < n; ++j) cin >> A[i][j];

        cout << "Entrez les valeurs du vecteur b (taille m) :" << endl;
        for(int i = 0; i < m; ++i) cin >> b[i];
            

            cout << "Entrez les valeurs du vecteur c (taille n) :" << endl;
        for(int i = 0; i < n; ++i) cin >> c[i];

        simplex<double>(A, b, c);
        return 0;
}
