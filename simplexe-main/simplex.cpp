#include <bits/stdc++.h>
using namespace std;

template <typename T>
void display(const vector<vector<T>> &A) {
        for(int i = 0; i < A.size(); ++i) {
                for(int j = 0; j < A[0].size(); ++j) cout << A[i][j] << ' ';
                cout << endl;
        }
        cout << endl;
}

template <typename T>
void display(const vector<T> &vec) {
        for(int i = 0; i < vec.size(); ++i) cout << vec[i] << ' ';
        cout << endl;
}

template<typename T>
void echange(vector<vector<T>> &mat, int i, int j) {
    for(int k = 0; k < mat.size(); ++k) {
        T tmp     = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = tmp;
    }
}

template<typename T>
void transpose(vector<vector<T>> &mat) {
    for(int i = 0;  i < mat.size(); ++i) {
        for(int j = 0; j < i; ++j) {
            auto tmp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = tmp;
        }
    }
}

template <typename T>
vector<vector<T>> inverse(vector<vector<T>> mat) {
    vector<vector<T>> ret(mat.size(), vector<T>(mat.size()));
    for(int i = 0; i < mat.size(); ++i) ret[i][i] = 1;
    bool pass = false;

        for(int i = 0; i < mat.size()-1; ++i) {
            //on cherche le pivot
            int p = -1;
            for(int j = i; j < mat.size(); ++j) {
                if(mat[i][j] != 0){
                    p = i;
                    break;
                }
            }
            if(p == -1) {//matrice non-inversible
                return vector<vector<T>>();
            }

            echange(mat, i, p);
            for(int k = i+1; k < mat.size(); ++k) {
                T alpha = mat[k][i] /mat[i][i];
                cout << alpha << endl;
                for(int l = 0; l < mat.size(); ++l) {
                    mat[k][l] -= alpha * mat[i][l];
                    ret[k][l] -= alpha * ret[i][l];
                }
                display(mat);
                cout << "end"<<endl;
            }
            display(ret);
        }


    display(ret);

}

template <typename T>
T euclidian_distance(vector<T> &v, vector<T> &u) {
        T ret = 0;
        for(int i = 0; i < v.size(); ++i) ret += pow((v[i] - u[i]), 2);
        return sqrt(ret);
}

template <typename T>
void pivot(vector<vector<T>> &A, vector<T> &b, vector<T> &c, T &nu, int e, int l) {
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
pair<int, int> pivot_choice_SE(const vector<vector<T>> &A, const vector<T> &b, const vector<T> &c, vector<int> &base, vector<int> &hbase) {
        //the pivot choice here is base on steepest edge
        //the pivot choice here is base on greatest improvement(GI)
        int entering     = -1;
        int outgoing     = -1;
        T eq_value       = numeric_limits<T>::min();

        vector<T> actual(A[0].size());
        for(int i = 0; i < A.size(); ++i) if(base[i] < A[0].size()) actual[base[i]] = b[i];
        //cout << "INITAIL : ";
        //display(actual);
        //cout << "CHOICE AMONG : " << endl;

        for(int i = 0; i < c.size(); ++i) {
                if(c[i] > 0) {
                        int eq = 0;
                        T eq_value_tmp = (A[0][i] < 0) ?  -(b[0]/A[0][i]) : numeric_limits<T>::max();
                        for(int k = 1; k < A.size(); ++k) {
                                if(A[k][i] >= 0) continue;
                                auto value = -(b[k] / A[k][i]);
                                if(value < eq_value_tmp) {
                                        eq_value_tmp       = value;
                                        eq                 = k;
                                }
                        }
                        if(eq_value_tmp == numeric_limits<T>::max()) {
                                cout <<"NON BORNE"<<endl;
                                return {-1,-1};
                        }

                        //cout << "WITH VALUE : " << eq_value_tmp << endl;
                        //cout << "VECTOR ... ";
                        vector<T> possibly_new(A[0].size());

                        auto tmp = -b[eq] / A[eq][i];
                        swap(hbase[i], base[eq]);
                        for(int k = 0; k < A.size(); ++k){
                                 if(base[k] < A[0].size()) {
                                         if(k == eq) possibly_new[base[k]] = tmp;
                                         else possibly_new[base[k]] = b[k] + tmp*A[k][i];//b[k]    += A[k][e] * b[l];
                                 }
                         }
                         //display(possibly_new);
                         swap(hbase[i], base[eq]);

                        auto step = eq_value_tmp/euclidian_distance(possibly_new, actual);
                        if(eq_value < step) {
                                eq_value = step;
                                outgoing = eq;
                                entering = i;
                        }
                }
        }
        return {entering, outgoing};
}

template <typename T>
pair<int, int> pivot_choice_GI(const vector<vector<T>> &A, const vector<T> &b, const vector<T> &c) {
        //the pivot choice here is base on greatest improvement(GI)
        int entering     = -1;
        int outgoing     = -1;
        T eq_value       = numeric_limits<T>::min();
        //cout << "CHOICE AMONG : " << endl;

        for(int i = 0; i < c.size(); ++i) {
                if(c[i] > 0) {
                        int eq = 0;
                        T eq_value_tmp = (A[0][i] < 0) ?  -(b[0]/A[0][i]) : numeric_limits<T>::max();
                        for(int k = 1; k < A.size(); ++k) {
                                if(A[k][i] >= 0) continue;
                                auto value = -(b[k] / A[k][i]);
                                if(value < eq_value_tmp) {
                                        eq_value_tmp       = value;
                                        eq                 = k;
                                }
                        }
                        if(eq_value_tmp == numeric_limits<T>::max()) {
                                cout <<"NON BORNE"<<endl;
                                return {-1,-1};
                        }

                        if(eq_value < eq_value_tmp) {
                                //cout << "WITH VALUE : " << eq_value_tmp << endl;
                                //cout << "VECTOR ... (" <<eq << ' ' << i << ')'<<endl;
                                eq_value = eq_value_tmp;
                                outgoing = eq;
                                entering = i;
                        }
                }
        }
        return {entering, outgoing};
}

template<typename T>
T dot(const vector<T> &w, const vector<T> &u, const vector<T> &v) {
    T ret = 0;
    for(int i = 0; i < w.size(); ++i) ret += w[i] * (u[i] - v[i]);
    return ret;
}

template <typename T>
pair<int, int> pivot_choice_MS(const vector<vector<T>> &A, const vector<T> &b, const vector<T> &c, const vector<T> &w, vector<int> &base, vector<int> &hbase) {
        //the pivot choice here is base on max slope
        int entering     = -1;
        int outgoing     = -1;
        T eq_value       = numeric_limits<T>::min();

        vector<T> v(A[0].size());
        for(int i = 0; i < A.size(); ++i) if(base[i] < A[0].size()) v[base[i]] = b[i];
    //    cout << "INITAIL : ";
    //    display(v);
    //    cout << "CHOICE AMONG : " << endl;

        for(int i = 0; i < c.size(); ++i) {
                if(c[i] > 0) {
                        int eq = 0;
                        T eq_value_tmp = (A[0][i] < 0) ?  -(b[0]/A[0][i]) : numeric_limits<T>::max();
                        for(int k = 1; k < A.size(); ++k) {
                                if(A[k][i] >= 0) continue;
                                auto value = -(b[k] / A[k][i]);
                                if(value < eq_value_tmp) {
                                        eq_value_tmp       = value;
                                        eq                 = k;
                                }
                        }
                        if(eq_value_tmp == numeric_limits<T>::max()) {
                                cout <<"NON BORNE"<<endl;
                                return {-1,-1};
                        }

                        vector<T> u(A[0].size());

                        auto tmp = -b[eq] / A[eq][i];
                        swap(hbase[i], base[eq]);
                        for(int k = 0; k < A.size(); ++k){
                                 if(base[k] < A[0].size()) {
                                         if(k == eq) u[base[k]] = tmp;
                                         else u[base[k]] = b[k] + tmp*A[k][i];//b[k]    += A[k][e] * b[l];
                                 }
                         }
                         swap(hbase[i], base[eq]);
                    //     cout << "DOT RESULT  = " <<  dot(w, u, v) << endl;
                    //     cout << "AND EQ_VALUE = " << eq_value_tmp << endl;

                        auto step = dot(w, u, v)/eq_value_tmp;
                    //    cout << "WITH VALUE : " << step << endl;
                    //    cout << "VECTOR ... ";
                    //    display(u);

                        if(eq_value < step) {
                                eq_value = step;
                                outgoing = eq;
                                entering = i;
                        }
                }
        }
        return {entering, outgoing};
}

template<typename T>
void initial_value(vector<T> &w, int n) {
    w.resize(n);
    for(int i = 0; i < n; ++i) w[i] = -1;
}


template <typename T>
vector<T> simplex(vector<vector<T>> A, vector<T> b, vector<T> c, int &step, vector<T> w = {}) {
        for(int i = 0; i < A.size(); ++i) for(int j = 0; j < A[0].size(); ++j) A[i][j] = -A[i][j];

        vector<int> hbase(A[0].size());
        vector<int> base(A.size());
        for(int i = 0; i < A[0].size(); ++i)    hbase[i]  = i;
        for(int i = 0; i < A.size(); ++i) base[i] = A[0].size()+i;

        if(w.empty()) initial_value(w, c.size());

        T nu = 0;
        pair<int, int> pivot_values = pivot_choice_MS(A,b,c, w,hbase, base);
        while(pivot_values.first >= 0) {
                //cout << endl;
                vector<T> optimal_tmp(A[0].size());
                for(int i = 0; i < A.size(); ++i) if(base[i] < A[0].size()) optimal_tmp[base[i]] = b[i];
                //display(optimal_tmp);
                //cout << "nu = " << nu << endl<<endl;
                //cout << "enter = " << pIvot_values.first << " out "<< pivot_values.second<<endl;
                //display current best solution with the associated vector

                pivot(A,b,c,nu,pivot_values.first,pivot_values.second);
               // display(A);
                swap(hbase[pivot_values.first], base[pivot_values.second]);
                pivot_values = pivot_choice_MS(A,b,c,w,hbase, base);
                ++step;
        }
        cout << "NUMBER OF STEP =  "<< step<<endl;
        vector<T> optimal(A[0].size());
        for(int i = 0; i < A.size(); ++i) if(base[i] < A[0].size()) optimal[base[i]] = b[i];

        display(optimal);
        cout << "nu = " << nu << endl;

        return optimal;
}

template <typename T, typename U, typename V>
struct triplet{
        T first;
        U second;
        V third;
};

#define N 1000
template<typename T>
vector<triplet<vector<T>, T, int>> test(const vector<vector<T>> &A, const vector<T> &b, const vector<T> &c) {
    vector<triplet<vector<T>, T, int>> ret(N);
    for(int i = 0; i < N; ++i) {
        vector<T> w = {cos((i*M_PI)/N), sin((i*M_PI)/N), i};
        int step = 0;
        auto max_vec = simplex(A,b,c,step,w);
       // cout << "\nFIN TEST : " << i+1 << endl<<endl;;

        T vl = 0;
        for(int i = 0; i < max_vec.size(); ++i) vl += max_vec[i]*c[i];
        //display(max_vec);
        ret[i] = {w, vl, step};
    }
    return ret;
}

template <typename T>
triplet<vector<vector<T>>, vector<T>, vector<T>> generate(int m, int n) {
    vector<vector<T>> a(m, vector<T>(n));
    vector<T> b(m);
    vector<T> c(n);
    
    for(int i = 0; i < m; ++i) {
        b[i] = rand();
        for(int j = 0; j < n; ++j) {
            a[i][j] = rand();
         //  cout << a[i][j] << ' ';
        }
        //cout << endl;
    }
    for(int j = 0; j < n; ++j) c[j] = rand();
    return {a, b, c};
}

int main() {
        //vector<vector<double>> A{{-1,-1},{-1,-2},{5,3}};
        //vector<double> c{1,1};90
        //vector<double> b{-5,-3,20};
        //vector<vector<double>> A{{15,20,25, 59, 57, 78},{35,60,60,35,60,60},{20,30,25,35,60,60},{0,250,0, 0,250,0}};
        //vector<double> c{300,250,450, 300, 250, 450};
        //vector<double> b{1200,3000,1500,500};
        vector<vector<double>> A{{1,0,0,},{20,1,0},{200,20,1}};//Klee minty example
        vector<double> c{100,10,1};
        vector<double> b{1,100,10000};
        //vector<vector<double>> A{{1,0},{20,1}};//Klee minty example
        //vector<double> c{100,10};
        //vector<double> b{1,100};
        int step = 0;
        //auto mx = simplex<double>(A,b,c, step);
        //display(mx);
    //    auto data = test(A,b,c);

        auto tr = generate<double>(10, 100);
        simplex<double>(tr.first, tr.second, tr.third, step);

        //for(int i = 0; i < data.size(); ++i) {
            //display((data[i]).first);
            //cout << "WITH VALUE : " << (data[i]).second<<endl;
          //  cout << "STEP NUMBER = " << (data[i]).third << endl;
        //}
        
        // shadow_vertex_method(A,b,c, x0, d, base);

       //vector<vector<double>> in{{1,2,3},{4,5,6},{7,8,9}};
      // display(in);
       //transpose(in);
       //display(in);
      // inverse(in);

        //vector<double> x0{1,1};
        //initialise(A, b, c, x0);
        //display(A);
        //display(c);
        return 0;
}
