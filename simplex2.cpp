#include <bits/stdc++.h>
#include "../../matrix/matrix.hpp"
using namespace std;

template <typename T>
void display(const vector<vector<T>> &A, const string& message = "") {
        if (!message.empty()) cout << message << endl;
        for(int i = 0; i < A.size(); ++i) {
                for(int j = 0; j < A[0].size(); ++j) cout << A[i][j] << ' ';
                cout << endl;
        }
        cout << endl;
}

template <typename T>
void display(const vector<T> &vec, const string& message = "") {
        if (!message.empty()) cout << message << endl;
        for(int i = 0; i < vec.size(); ++i) cout << vec[i] << ' ';
        cout << endl;
}

void enter(vector<vector<double>> &A, vector<double>&b, vector<double> &c) {
        cout << "Veuillez entrer les coefficients de la matrice A (contraintes) :" << endl;
        for(int i = 0; i < A.size(); ++i) 
            for(int j = 0; j < A[i].size(); ++j) 
                cin >> A[i][j];
                
        cout << "Veuillez entrer le vecteur b (valeurs constantes des contraintes) :" << endl;
        for(int i = 0; i < b.size(); ++i) 
            cin >> b[i];
            
        cout << "Veuillez entrer le vecteur c (coefficients de la fonction objectif) :" << endl;
        for(int i = 0; i < c.size(); ++i) 
            cin >> c[i];
}

void compute_current_position(vector<double> &position, const vector<vector<double>> &A, const vector<double> &b, const vector<int> &base) {
        for(int i = 0; i < position.size(); ++i) 
            position[i] = 0;

        Matrix<double> B(base.size(), base.size());
        for(int i = 0; i < base.size(); ++i) 
            for(int j = 0; j < base.size(); ++j) 
                B(i, j) = A[i][base[j]];

        cout << "Affichage de la matrice B (base actuelle) :" << endl;
        B.display();

        auto inv = B.gaussian_inverse();
        cout << "Affichage de l'inverse de B :" << endl;
        inv.display();

        display(b, "Affichage du vecteur b :");
        auto tmp = inv.dot(b);

        for(int i = 0; i < base.size(); ++i)
            position[base[i]] = tmp[i];
}

bool continue_simplex(const vector<double> &c, const vector<int> &base) {
        for(int i = 0; i < base.size(); ++i) 
            if (c[base[i]] > 0) 
                return true;
        return false;
}

void simplex(vector<vector<double>> &A, vector<double> &b, vector<double> &c) {
        int m = A.size();
        int n = c.size();
        vector<int> base(m);
        vector<int> hbase(n);

        for(int i = 0; i < m; ++i) {
                A[i].resize(A[i].size() + m);
                A[i][m + i] = 1;
                base[i] = m + i;
        }

        for(int i = 0; i < n; ++i) 
            hbase[i] = i;

        for(int i = 0; i < m; ++i) 
            c.push_back(0);

        vector<double> current_position(n + m);
        
        display(A, "Affichage de la matrice augmentée A :");
        compute_current_position(current_position, A, b, base);
        
        display(current_position, "Affichage de la position actuelle :");

        while (continue_simplex(c, base)) {
                for (int i = 0; i < base.size(); ++i) {
                        if (c[base[i]] < 0) continue;
                        double mn = numeric_limits<double>::max();
                        int ele_swap = -1;

                        for (int j = 0; j < hbase.size(); ++j) {
                                if (A[base[i]][hbase[j]] >= 0) continue;
                                auto value = -(b[base[i]] / A[base[i]][hbase[j]]);
                                if (value < mn) {
                                        mn = value;
                                        ele_swap = i;
                                }
                        }
                }
        }

        // Affichage de la solution optimale
        cout << "LA solution optimale est : " << endl;
        display(current_position);
}

int main() {
        int n, m; // m le nombre de contraintes et n le nombre de variables
        cout << "Veuillez entrer le nombre de contraintes (m) et de variables (n) :" << endl;
        cin >> m >> n;

        vector<vector<double>> A(m, vector<double>(n));
        vector<double> b(m);
        vector<double> c(n);

        enter(A, b, c);
        simplex(A, b, c);

        return 0;
}

