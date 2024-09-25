#include <iostream>
#include <unistd.h>
#include <tuple>
#include <limits>
#include <vector>
#include <cmath>
#include <list>
#include <queue>
#include <algorithm>
#include <set>
#include <cmath>
#include <functional>
#include <fstream>
#include <memory>
#include <map>
#include <assert.h>

template <typename T> class List;
template <typename T> class Node;
template <typename T> class Node_Closure;
template <typename T> class DT;


namespace Util {
        template <typename T>
        T min(const T& a, const T& b) {
                return (a < b) ? a : b;
        }

        template <typename T>
        T sum(const std::vector<T> &v) {
                int ret = 0;
                for(int i = 0; i < v.size(); i++)
                        ret += v[i];
                return ret;
        }

        template <typename T>
        T max(const T& a, const T& b) {
                return (a > b) ? a : b;
        }

        template <typename T>
        void swap_(T &a, T &b) {
                T tmp = a;
                a = b;
                b = tmp;
        }

        template <typename T>
        void display(T *a, int s) {
                for(int i = 0; i < s; i++)
                        std::cout << a[i] << " ";
                std::cout << std::endl;
        }

        template <typename T>
        void display(const std::vector<T> &a) {
                for(int i = 0; i < a.size(); i++)
                        std::cout << a[i] << " ";
                std::cout << std::endl;
        }

        template <typename T>
        void display(const std::vector<std::set<T>> &v) {
                for(int i = 0; i < v.size(); i++) {
                        for(auto &e : v[i])
                                std::cout << e << " ";
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<std::vector<T>> &a) {
                for(int i = 0; i < a.size(); i++)
                        Util::display(a[i]);
        }

        template <typename T>
        void generate(T *tab, int n) {
                srand(getpid());
                for(int i = 0; i < n; i++)
                        tab[i] = rand() % 100;
        }

        template <typename T>
        void generate(std::vector<T> &tab) {
                srand(getpid());
                for(int i = 0; i < tab.size(); i++)
                        tab[i] = rand() % 100;
        }

        template <typename T>
        T min(T *tab, int n) {
                T min = tab[0];
                for(int i = 1; i < n; i++) {
                        if(tab[i] < min)
                                min = tab[i];
                }
                return min;
        }

        template <typename T>
        T max(T *tab, int n) {
                T max = tab[0];
                for(int i = 1; i < n; i++) {
                        if(tab[i] > max)
                                max = tab[i];
                }
                return max;
        }

        template <typename T>
        std::pair<T, T> max_min(T *a, int n) {
                std::pair<T, T> ret(a[0], a[0]);
                for(int i = 1; i < n; i+=2) {
                        int max = i;
                        int min = i;
                        if(a[i] > a[i-1])
                                min = i-1;
                        if(a[max] > ret.first)
                                ret.first = a[max];
                        if(a[min] < ret.second)
                                ret.second = a[min];
                }
                if(n % 2 == 1) {
                        if(a[n-1] > ret.first)
                                ret.first = a[n-1];
                        else if(a[n-1] < ret.second)
                                ret.second = a[n-1];
                }
                return ret;
        }

        template <typename T, typename U>
        bool pair_comparaison_dikjstra(const std::pair<T, U> &pr1, const std::pair<T, U> &pr2) {//we compare the second element
                return pr1.second < pr2.second;
        }

        template <typename T>
        void display2D(const std::vector<T> &vec) {
                int n = sqrt(vec.size());
                for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++)
                                std::cout << vec[i*n+j] << " ";
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<std::list<T>> &gr) {
                for(int i = 0; i < gr.size(); i++) {
                        std::cout << "Vertices " << i + 1 << " : " ;
                        for(auto e : gr[i])
                                std::cout << e << " ";
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<std::pair<T, T>> &vec) {
                for(int i = 0; i < vec.size(); i++) {
                        std::cout << "NODE " << i << " : " << vec[i].first << " and " << vec[i].second << std::endl;
                }
        }

        void display(const std::vector<std::list<std::pair<int, int>>> &gr) {
                for(int i = 0; i < gr.size(); i++) {
                        std::cout << "Vertices : " << i + 1 << std::endl;
                        for(auto e : gr[i])
                                std::cout << e.first << " value : " << e.second << std::endl;
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(const std::vector<List<std::pair<T, T>>> &v) {
                for(int i = 0; i < v.size(); i++) {
                        std::cout << "Vertices : " << i + 1 << std::endl;
                        auto it = v[i].get();
                        for(int j = 0; j < v[i].size(); j++) {
                                if(it == nullptr)
                                        break;
                                std::cout << it->value().first << " value : " << it->value().second << std::endl;
                                it = it->next();
                        }
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(std::vector<List<T>> &v) {
                for(int i = 0; i < v.size(); i++) {
                        std::cout << "Vertices : " << i + 1 << std::endl;
                        auto it = v[i].get();
                        for(int j = 0; j < v[i].size(); j++) {
                                if(it == nullptr)
                                        break;
                                std::cout << it->value() << " " << std::endl;
                                it = it->next();
                        }
                        std::cout << std::endl;
                }
        }

        template <typename T>
        void display(List<T> &ls) {
                auto it = ls.get();
                for(int i = 0; i < ls.size(); i++) {
                        std::cout << it->value() << std::endl;
                        it = it->next();
                }
        }

        template <typename T>
        void display(Node_Closure<T> &ls) {
                if(ls.is_list())
                        ls.get_list()->display();
                else if(ls.is_node()) {
                        std::cout << ls.get_node()->value() << " AND SIZE = " << ls.size() << std::endl;
                }
        }

        template <typename T>
        void display(std::vector<Node_Closure<T>> &v) {
                for(int i = 0; i < v.size(); i++) {
                        Util::display(v[i]);
                }
        }

        template <typename T>
        void display(std::vector<List<DT<T>>> &v) {
                for(int i = 0; i < v.size(); i++) {
                        auto it = v[i].get();
                        for(int j = 0; j < v[i].size(); j++) {
                                if(it == nullptr)
                                        break;
                                it->value().display();
                                it= it->next();
                        }
                }
        }

        struct triplet{
                int x;
                int y;
                int z;
        };
        typedef struct triplet triplet;

        bool is_in(std::vector<std::set<int>> &en, int i, int s) {
                for(auto e : en[i]) {
                        if(e == s)
                                return true;
                }
                return false;
        }

        template <typename T>
        bool is_in(std::list<T> &ls, T ele) {
                for(auto &e : ls)
                        if(e == ele)
                                return true;
                return false;
        }

        template <typename T>
        bool is_in(std::set<T> &en, T ele) {
                return (en.find(ele) == en.end()) ? false : true;
        }

        template <typename T>
        bool is_in(List<T> &ls, T ele) {
                auto it  = ls.get();
                for(int i = 0; i < ls.size(); i++) {
                        if(it == nullptr)
                                return false;
                        if(it->value() == ele)
                                return true;
                        it = it->next();
                }
                return false;
        }

        template <typename T>
        bool equal(const std::vector<std::vector<T>> &gr1, const std::vector<std::vector<T>> &gr2) {
                for(int i = 0; i < gr1.size(); i++) {
                        for(int j = 0; j < gr1.size(); j++) {
                                if(gr1[i][j] != gr2[i][j])
                                        return false;
                        }
                }
                return true;
        }

};

namespace Matrix {
        template <typename T>
        std::vector<std::vector<T>> dot(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {//naive multiplication O(n^3)
                if(a[0].size() != b.size())
                        throw std::invalid_argument("bad size for matrix multiplication\n");

                std::vector<std::vector<T>> ret(a.size(), std::vector<T>(b[0].size()));
                std::vector<T> tmp(b.size());

                for(int i = 0; i < b[0].size(); i++) {
                        for(int j = 0; j < a[0].size(); j++)//cache-friendly
                                tmp[j] = b[j][i];

                        for(int j = 0; j < a.size(); j++) {
                                T sum = 0;
                                for(int k = 0; k < a[0].size(); k++) {
                                        sum += a[j][k] * tmp[k];
                                }
                                ret[j][i] = sum;
                        }
                }
                return ret;
        }

        template <typename T>
        std::vector<std::vector<T>> add(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {//O(n^2)
                std::vector<std::vector<T>> ret(a.size(), std::vector<T>(a[0].size()));
                for(int i = 0; i < a.size(); i++) {
                        for(int j = 0; j < a[0].size(); j++) {
                                ret[i][j] = a[i][j] + b[i][j];
                        }
                }
                return ret;
        }

        template <typename T>
        std::vector<std::vector<T>> sub(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {//O(n^2)
                std::vector<std::vector<T>> ret(a.size(), std::vector<T>(a[0].size()));
                for(int i = 0; i < a.size(); i++) {
                        for(int j = 0; j < a[0].size(); j++) {
                                ret[i][j] = a[i][j] - b[i][j];
                        }
                }
                return ret;
        }

        template <typename T>
        std::vector<std::vector<T>> strassen(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {//O(n^lg2(7))
                if(a[0].size() != b.size())
                        throw std::invalid_argument("bad size for matrix multiplication\n");
                if(a.size() <= 1)
                        return Matrix::dot(a, b);
                auto n = a.size() >> 1;
                std::vector<std::vector<T>> a11(n, std::vector<T>(n));
                std::vector<std::vector<T>> a12(n, std::vector<T>(n));
                std::vector<std::vector<T>> a21(n, std::vector<T>(n));
                std::vector<std::vector<T>> a22(n, std::vector<T>(n));

                std::vector<std::vector<T>> b11(n, std::vector<T>(n));
                std::vector<std::vector<T>> b12(n, std::vector<T>(n));
                std::vector<std::vector<T>> b21(n, std::vector<T>(n));
                std::vector<std::vector<T>> b22(n, std::vector<T>(n));

                for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++) {
                                a11[i][j] = a[i][j];
                                a12[i][j] = a[i + n][j];
                                a21[i][j] = a[i][j + n];
                                a22[i][j] = a[i + n][j + n];

                                b11[i][j] = b[i][j];
                                b12[i][j] = b[i + n][j];
                                b21[i][j] = b[i][j + n];
                                b22[i][j] = b[i + n][j + n];
                        }
                }

                auto m1  = strassen(add(a11 , a22) , add(b11, b22));
                auto m2  = strassen(add(a21, a22) , b11);
                auto m3  = strassen(a11, sub(b12, b22));
                auto m4  = strassen(a22, sub(b21, b11));
                auto m5  = strassen(add(a11, a12) , b22);
                auto m6  = strassen(sub(a21, a11) , add(b11, b12));
                auto m7  = strassen(sub(a12, a22) , add(b21, b22));

                auto c11  = add(sub(add(m1, m4), m5), m7);
                auto c12  = add(m3, m5);
                auto c21  = add(m2, m4);
                auto c22  = add(add(sub(m1, m2), m3), m6);

                std::vector<std::vector<T>> ret(a.size(), std::vector<T>(a.size()));
                for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++) {
                                ret[i][j] = c11[i][j];
                                ret[i + n][j] = c12[i][j];
                                ret[i][j + n] = c21[i][j];
                                ret[i + n][j + n] = c22[i][j];
                        }
                }
                return ret;
        }
};


#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) + 1)
#define PARENT(i) (i >> 1)

template <typename T>
class Heap {//max
public :
        Heap(const std::vector<T> &vec) : tab_(vec.size()) {
                std::copy(vec.begin(), vec.end(), this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        Heap(T *tab, int s) : tab_(s){
                std::copy(tab, tab + s, this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        void insert(T ele, int p) {//to insert element at the end and get it up
                if(p >= this->tab_.size())
                        this->tab_.push_back(p);
                else
                        this->tab_[p] = ele;

                auto p2 = PARENT(p);
                while(this->tab_[p] > this->tab_[p2] && p2 < p) {
                        Util::swap_(this->tab_[p], this->tab_[p2]);
                        p = p2;
                        p2 = PARENT(p);
                }
        }

        void insert(int i) {//insert the element at index i
                if(i >= this->tab_.size() || i < 0)
                        return ;
                int l   = LEFT(i);
                int r   = RIGHT(i);
                int max = i;
                if(l < this->tab_.size() && this->tab_[i] < this->tab_[l])
                        max = l;
                if(r < this->tab_.size() && this->tab_[max] < this->tab_[r])
                        max = r;
                if(i == max)
                        return ;

                Util::swap_(this->tab_[i], this->tab_[max]);
                this->insert(max);
        }

        T pop() {
                T ret         = this->tab_[0];
                this->tab_[0] = this->tab_.back();
                this->tab_.pop_back();
                this->insert(0);
                return ret;
        }

        T max() {
                return this->tab_[0];
        }

        void display() {
                Util::display(this->tab_);
        }

private :
        std::vector<T> tab_;
};

template <typename T>
class Heap_min {
public :
        Heap_min(const std::vector<T> &vec) : tab_(vec.size()) {
                std::copy(vec.begin(), vec.end(), this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        Heap_min(T *tab, int s) : tab_(s){
                std::copy(tab, tab + s, this->tab_.begin());
                for(int i = (this->tab_.size() >> 1); i >= 0; i--) {
                        this->insert(i);
                }
        }

        void insert(T ele, int p) {//to insert element at the end and get it up
                if(p >= this->tab_.size())
                        this->tab_.push_back(p);
                else
                        this->tab_[p] = ele;

                auto p2 = PARENT(p);
                while(this->tab_[p] < this->tab_[p2] && p2 > p) {
                        Util::swap_(this->tab_[p], this->tab_[p2]);
                        p  = p2;
                        p2 = PARENT(p);
                }
        }

        void insert(int i) {//insert the element at index i
                if(i >= this->tab_.size() || i < 0)
                        return ;
                int l   = LEFT(i);
                int r   = RIGHT(i);
                int max = i;
                if(l < this->tab_.size() && this->tab_[i] > this->tab_[l])
                        max = l;
                if(r < this->tab_.size() && this->tab_[max] > this->tab_[r])
                        max = r;
                if(i == max)
                        return ;

                Util::swap_(this->tab_[i], this->tab_[max]);
                this->insert(max);
        }

        T pop() {
                T ret         = this->tab_[0];
                this->tab_[0] = this->tab_.back();
                this->tab_.pop_back();
                this->insert(0);
                return ret;
        }

        T min() {
                return this->tab_[0];
        }

        void display() {
                Util::display(this->tab_);
        }

private :
        std::vector<T> tab_;
};

template <typename T> class Fibonacci_heap;

template <typename T>
class Fibonacci_root {
  public:
        Fibonacci_root(T value = 0, int n = 0, Fibonacci_root<T> *p = nullptr, Fibonacci_root<T> *anild = nullptr, Fibonacci_root<T> *left = nullptr, Fibonacci_root<T> *right = nullptr) : value_(value),  n_(n){
                this->p_      = p;
                this->child_  = child_;
                this->left_   = left;
                this->right_  = right;
                this->marked_ = false;
        }

        Fibonacci_root(const Fibonacci_root<T>& nd) : value_(nd.value_), n_(nd.n_) {
                this->p_       = nd.p_;
                this->child_   = nd.child_;
                this->left_    = nd.left_;
                this->right_   = nd.right_;
                this->marked_  = nd.marked_;
        }

        T get_value() {
                return this->value_;
        }

        ~Fibonacci_root() {
                if(this->p_ != nullptr)
                        delete this->p_;
                if(this->child_ != nullptr)
                        delete this->child_;
                if(this->left_ != nullptr)
                        delete this->left_;
                if(this->right_ != nullptr)
                        delete this->right_;
        }

        void display() {
                auto tmp = this;
                std::cout << "VALUE = " << tmp->value_ << std::endl;
                auto it  = tmp->right_;
                if(this->child_ != nullptr) {
                        std::cout << "SON : " <<this->child_<< std::endl;
                        this->child_->display();
                        std::cout << "ON REMONTE DUN NIVEAU" << std::endl;
                }else
                        std::cout << "SON  = NULLPTR" << std::endl;

                while(it != tmp) {
                        std::cout << "VALUE = " << it->value_ << std::endl;
                        if(it->child_ != nullptr) {
                                std::cout << "SON : "  << std::endl;
                                it->child_->display();
                                std::cout << "ON REMONTE DUN NIVEAU" << std::endl;
                        }else
                                std::cout << "SON = NULLPTR" << std::endl;

                        it = it->right_;
                }
        }

        void display_dikjstra() {
                auto it = this;
                std::cout << "FIRST = " << it->value_.first << std::endl;
                std::cout << "SECOND = " << it->value_.second << std::endl;
                if(this->child_ != nullptr) {
                        std::cout << "SON : " << std::endl;
                        this->child_->display_dikjstra();
                        std::cout << "ON REMONTE D UN NIVEAU" << std::endl;
                }else
                        std::cout << "SON = NULLPTR" << std::endl;
                it = it->right_;

                while(it != this) {
                        std::cout << "FIRST = " << it->value_.first << std::endl;
                        std::cout << "SECOND = " << it->value_.second << std::endl;
                        if(it->child_ != nullptr) {
                                std::cout << "SON : " << std::endl;
                                it->child_->display_dikjstra();
                                std::cout << "ON REMONTE D UN NIVEAU" << std::endl;
                        }else
                                std::cout << "SON = NULLPTR" << std::endl;
                        it = it->right_;
                }
        }

  private:
        T value_;
        int n_;
        Fibonacci_root<T> *p_;
        Fibonacci_root<T> *child_;
        Fibonacci_root<T> *right_;
        Fibonacci_root<T> *left_;
        bool marked_;
        friend class Fibonacci_heap<T>;
};

template <typename T>
class Fibonacci_heap {
  public:
          Fibonacci_heap(std::function<bool(T&, T&)> fct) : min_(nullptr), n_(0), function_(fct) {
          }

          Fibonacci_root<T>* insert(T ele) {
                auto n = new Fibonacci_root<T>(ele);
                this->insert_root(n);
                return n;
          }

          void insert_root(Fibonacci_root<T> *n) {
                  if(this->min_ == nullptr) {
                          this->min_         = n;
                          this->min_->right_ = n;
                          this->min_->left_  = n;
                  }
                  else{
                          n->left_           = this->min_;
                          n->right_          = this->min_->right_;
                          this->min_->right_ = n;
                          n->right_->left_   = n;

                          if(this->function_(n->value_, this->min_->value_))
                                  this->min_ = n;
                  }
                  ++this->n_;
          }

          void remove_root(Fibonacci_root<T> *&n) {//don t care of child_ attribute
                if(n->p_ != nullptr)
                        --n->p_->n_;
                else
                        --this->n_;
                n->child_        = nullptr;
                n->right_->left_ = n->left_;
                n->left_->right_ = n->right_;
                n->left_         = nullptr;
                n->right_        = nullptr;
                n->p_            = nullptr;
          }

          void union_heap(Fibonacci_heap<T> &hp) {
                  if(hp.min_ == nullptr)
                        return ;
                  auto tmp = hp.min_;

                  tmp->left_->right_              = this->min_->right_;
                  this->min_->right_->left_       = tmp->left_;
                  tmp->left_                      = this->min_;
                  this->min_->right_              = tmp;

                  if(this->min_ == nullptr || (hp.min_ != nullptr && this->function_(hp.min_->value_, this->min_->value_)))
                        this->min_ = hp.min_;
                this->n_ += hp.n_;
          }

          void link(Fibonacci_root<T> *a, Fibonacci_root<T> *b) {
                a->marked_               = false;
                b->left_->right_         = b->right_;
                b->right_->left_         = b->left_;

                if(a->child_ == nullptr) {
                        a->child_ = b;
                        b->left_  = b;
                        b->right_ = b;
                }
                else {
                        b->right_                = a->child_->right_;
                        a->child_->right_->left_ = b;
                        b->left_                 = a->child_;
                        a->child_->right_        = b;
                }
                b->p_ = a;
                --this->n_;
                ++a->n_;
          }

          bool empty() {
                  return (this->min_ == nullptr) ? true : false;
          }

          void cascade_cut(Fibonacci_root<T> *rt) {
                  auto p = rt->p_;
                  if(p != nullptr) {
                        if(!rt->marked_)
                                rt->marked_ = true;
                        else{
                                this->cut(rt);
                                this->cascade_cut(p);
                        }
                  }
          }

          void cut(Fibonacci_root<T> *rt) {
                auto p = rt->p_;
                if(p == nullptr)
                          return ;

                if(rt->left_ != rt) {

                        rt->left_->right_ = rt->right_;
                        rt->right_->left_ = rt->left_;
                }
                p->child_ = nullptr;
                this->insert_root(rt);
                --p->n_;
                rt->p_ = nullptr;
          }

          Fibonacci_root<T>* get_min() {
                  return this->min_;
          }

          void reduce(Fibonacci_root<T> *rt, T new_k) {
                if(this->function_(rt->value_, new_k))
                        throw std::invalid_argument("NEW KEY TOO BIG : \n");
                auto p     = rt->p_;
                rt->value_ = new_k;

                if(p != nullptr && this->function_(rt->value_, p->value_)) {
                        this->cut(rt);
                        this->cascade_cut(p);
                }
                if(this->function_(rt->value_, this->min_->value_))
                        this->min_ = rt;
          }

          void consolidate() {
                  auto tmp_n = this->n_;
                  int s = 0;
                  while (tmp_n >>= 1) ++s;
                  ++s;
                  ++s;

                  std::vector<Fibonacci_root<T>*> a(s, nullptr);
                  auto it                = this->min_->right_;
                  a[this->min_->n_]      = this->min_;

                  while(it != this->min_) {//premier
                          auto d      = it->n_;
                          auto x      = it;
                          auto future = it->right_;

                          while(a[d] != nullptr && a[d] != it) {
                                Fibonacci_root<T> *y = a[d];
                                if(this->function_(y->value_, x->value_)) {
                                        auto sw = x;
                                        x       = y;
                                        y       = sw;
                                }//x has the smallest key
                                this->link(x, y);
                                if(this->min_ == y)
                                        this->min_ = x;

                                a[d] = nullptr;
                                if(d < s-1)
                                        ++d;
                                else
                                        break;
                          }
                          it           = x;
                          a[d]         = it;
                          it           = it->right_;
                  }

                  auto ite = this->min_;
                  for(int i = 0; i < this->n_; i++) {
                        if(this->function_(ite->value_, this->min_->value_))
                                this->min_ = ite;
                        ite = ite->right_;
                  }
          }

          T pop() {//return min
                if(this->min_ == nullptr)
                        throw std::invalid_argument("CANT POP, HEAP = NULLPTR\n");

                auto ret         = this->min_->value_;
                if(this->min_->child_ != nullptr) {
                        auto itt = this->min_->child_->right_;
                        this->min_->child_->p_ = nullptr;

                        auto it = itt;
                        this->insert_root(this->min_->child_);
                        while(it != this->min_->child_) {

                                auto future = it->right_;
                                this->insert_root(it);
                                it->p_ = nullptr;
                                it     = future;
                                ++this->n_;
                        }
                        this->min_->child_ = nullptr;
                }

                if(this->min_ == this->min_->right_) {
                          this->min_->right_ = nullptr;
                          this->min_->left_  = nullptr;
                          this->min_->p_     = nullptr;
                          this->min_->child_ = nullptr;
                          delete this->min_;
                          this->min_         = nullptr;
                          this->n_           = 0;
                          return ret;
                }

                auto new_min = this->min_->right_;

                this->min_->child_        = nullptr;
                this->min_->right_->left_ = this->min_->left_;
                this->min_->left_->right_ = this->min_->right_;
                this->min_->left_         = nullptr;
                this->min_->right_        = nullptr;
                this->min_->p_            = nullptr;
                auto tmp_min              = this->min_;

                this->min_   = new_min;
                delete tmp_min;

                this->consolidate();
                return ret;
          }

          void debug_pointer() {
                if(this->min_ == nullptr)
                        return ;
                auto it = this->min_;
                std::cout << "\nBEGIN DEBUG POINTER" << std::endl;
                std::cout << it->value_ << std::endl;
                std::cout << (long int) it << std::endl;
                std::cout << (long int) it->left_ << std::endl;
                std::cout << (long int) it->right_ << std::endl;
                it = it->right_;

                while(it != this->min_) {
                        std::cout << "NEXT\n";
                        std::cout << it->value_ << std::endl;
                        std::cout << (long int) it << std::endl;
                        std::cout << (long int) it->left_ << std::endl;
                        std::cout << (long int) it->right_ << std::endl;
                        it = it->right_;
                }
                std::cout << "END DEBUG POINTER\n" << std::endl;
          }

          void display_dikjstra() {
                if(this->min_ == nullptr)
                        return ;
                this->min_->display_dikjstra();
          }

          void display() {
                if(this->min_ != nullptr)
                        this->min_->display();
                else
                        std::cout << "HEAP = NULLPTR" << std::endl;
          }

  private:
        Fibonacci_root<T> *min_;
        int n_;
        std::function<bool(T&, T&)> function_;
};

template <typename T> class BinarySearchTree;

template <typename T>
class node_tree {
public:
        node_tree(T value = 0) : key_(value), left_(nullptr), right_(nullptr) {
        }

        ~node_tree() {
                if(this->left_ != nullptr)
                        delete this->left_;
                if(this->right_ != nullptr)
                        delete this->right_;
        }

        T key() {
                return this->key_;
        }

        void add_left(T ele) {
                left_ = new node_tree<T>(ele);
        }

        node_tree<T> *left() {
                return this->left_;
        }

        node_tree<T> *right() {
                return this->right_;
        }

        void add_right(T ele) {
                right_ = new node_tree<T>(ele);
        }

        node_tree<T> *cut_right() {
                auto tmp = right_;
                right_   = nullptr;
                return tmp;
        }

        node_tree<T> *cut_left() {
                auto tmp = right_;
                right_   = nullptr;
                return tmp;
        }

        void set_key(T ele) {
                key_ = ele;
        }

private:
        T key_;
        node_tree<T> *left_;
        node_tree<T> *right_;
        friend class BinarySearchTree<T>;

};

template <typename T>
class BinarySearchTree {
public:
        BinarySearchTree(const std::vector<T> &vec) : nd_(nullptr){
                for(int i = 0; i < vec.size(); i++) {
                        this->insert(vec[i]);
                }
        }

        ~BinarySearchTree() {
                if(!nd_)
                        delete this->nd_;
        }

        T max() {
                auto tmp = this->nd_;
                while(tmp->right_ != nullptr)
                        tmp = tmp->right_;
                return tmp->key_;
        }

        T min() {
                auto tmp = this->nd_;
                while(tmp->left_ != nullptr)
                        tmp = tmp->left_;
                return tmp->key_;
        }

        node_tree<T> *search(T ele) {
                auto tmp = this->nd_;
                while(tmp->key_ != ele && tmp != nullptr) {
                        if(ele > tmp->key_)
                                tmp = tmp->right_;
                        else
                                tmp = tmp->left_;
                }
                return tmp;
        }

        void insert(T ele) {
                if(this->nd_ == nullptr) {
                        this->nd_ = new node_tree<T>(ele);
                        return ;
                }

                 auto tmp = this->nd_;
                 node_tree<T> *last;
                 while(tmp != nullptr) {
                         last = tmp;
                         if(ele > tmp->key_)
                                tmp = tmp->right_;
                        else if(ele < tmp->key_)
                                tmp = tmp->left_;
                        else
                                return ;
                 }
                 if(ele > last->key_) {
                        last->right_ = new node_tree<T>(ele);
                }else {
                        last->left_ = new node_tree<T>(ele);
                }
        }

private:
        node_tree<T> *nd_;
};

template <typename T>
class DT {
private:
        T vl_;
        std::pair<std::shared_ptr<Node<T>>, int> nd_;
public :
        DT(T value = 0) : vl_(value), nd_(std::pair<std::shared_ptr<Node<T>>, int>(nullptr, 0)){}
        bool is_value() {
                return (nd_.first == nullptr) ? true : false;
        }

        bool is_node() {
                return !is_value();
        }

        void set_value(T ele) {
                vl_ = ele;
        }

        T get_value() {
                if(!is_value())
                        throw std::invalid_argument("ERROR IT IS NOT A VALUE\n");
                return vl_;
        }

        void set_node(std::shared_ptr<Node<T>>nd, int sz) {
                nd_ = std::pair<std::shared_ptr<Node<T>>, int>(nd, sz);
        }

        void display() {
                if(is_value())
                        std::cout << vl_ << std::endl;
                else
                        std::cout << "IT S A POINTER TO : " << nd_.first.get() << " AND SIZE = " << nd_.second << std::endl;
        }

};

template <typename T>
class Node {
  private:
        std::shared_ptr<Node<T>> left_;
        std::shared_ptr<Node<T>> right_;
        T value_;
        friend class List<T>;
  public :
        Node(T value = 0) : left_(nullptr), right_(nullptr), value_(value) {}

        ~Node() {}

        T& value() {
                return value_;
        }

        std::shared_ptr<Node> next() {
                return right_;
        }

};

template <typename T>
class List {
  private:
        std::shared_ptr<Node<T>> begin_;
        std::shared_ptr<Node<T>> end_;
        int size_;
  public:
        List() : begin_(nullptr), end_(nullptr), size_(0) {}

        ~List() {}

        void push_back(T value) {
                if(size_ == 0) {
                        ++size_;
                        begin_ = std::make_shared<Node<T>>(value);
                        return ;
                }
                if(size_ == 1) {
                        ++size_;
                        end_                 = std::make_shared<Node<T>>(value);
                        begin_->right_       = end_;
                        end_->left_          = begin_;
                        return ;
                }

                ++size_;
                auto bc     = end_;
                end_        = std::make_shared<Node<T>> (value);
                bc->right_  = end_;
                end_->left_ = bc;
        }

        void concatenate(List &nd) {//at the end
                if(nd.size_ == 0)
                        return;

                if(size_ == 1) {
                        auto lt                = nd.begin_->left_;
                        size_                 += nd.size_;
                        begin_->right_         = nd.begin_;
                        nd.begin_->left_       = begin_;
                        if(lt != nullptr) {
                                begin_->left_  = lt;
                                lt->right_     = begin_;
                        }
                        if(nd.size_ == 1)
                                end_ = nd.begin_;
                        else
                                end_ = nd.end_;

                        return ;
                }
                if(nd.size_ == 1) {
                        ++size_;
                        nd.begin_->left_ = end_;
                        end_->right_ = nd.begin_;
                        return ;
                }
                size_                 += nd.size_;
                nd.begin_->left_       = end_;
                end_->right_           = nd.begin_;

                end_                   = nd.end_;
        }

        std::shared_ptr<Node<T>> get() {
                return begin_;
        }

        void display(){
                Util::display(*this);
        }

        int size() {
                return size_;
        }

        void push_front(T value) {
                if(size_ == 0) {
                        ++size_;
                        begin_ = std::make_shared<Node<T>>(value);
                        return ;
                }
                if(size_ == 1) {
                        ++size_;
                        end_              = std::make_shared<Node<T>>(value);
                        end_->left_       = begin_;
                        begin_->right_    = end_;
                }
                ++size_;
                auto tmp       = begin_;
                begin_         = std::make_shared<Node<T>>(value);
                begin_->right_ = tmp;
                tmp->left_     = begin_;
        }

};

template <typename T>
class Node_Closure {
private :
        std::shared_ptr<Node<T>> nd_;
        std::shared_ptr<List<T>> ls_;
        int sz_;
public :
        Node_Closure() : nd_(nullptr), ls_(nullptr), sz_(0) {}

        void set_list(std::shared_ptr<List<T>> &ls) {
                if(nd_ != nullptr)
                        throw std::invalid_argument("IS IT REALLY AN ACYCLIC GRAPH??\n");
                ls_ = ls;
        }

        std::shared_ptr<List<T>> retire_list() {
                if(ls_ == nullptr)
                        throw std::invalid_argument("AIE AIE AIE!!!\n");
                nd_         = ls_->get();
                auto ret    = ls_;
                sz_         = ls_->size();
                ls_         = nullptr;
                return ret;
        }

        void push_front(T ele) {
                if(ls_ == nullptr && nd_ == nullptr) {
                        ls_ = std::make_shared<List<T>>();
                        ls_->push_back(ele);
                        return ;
                }
                if(ls_ == nullptr)
                        return ;
                ls_->push_front(ele);
        }

        void push_back(T ele) {
                if(ls_ == nullptr && nd_ == nullptr) {
                        ls_ = std::make_shared<List<T>>();
                        ls_->push_back(ele);
                        return ;
                }
                if(ls_ == nullptr)
                        return ;
                ls_->push_back(ele);
        }

        bool is_list() {
                return (ls_ == nullptr) ? false : true;
        }

        std::shared_ptr<List<T>> get_list() {
                return ls_;
        }

        std::shared_ptr<Node<T>> get_node() {
                return nd_;
        }

        bool is_node() {
                return (nd_ == nullptr) ? false : true;
        }

        int size() {
                return sz_;
        }

        void add_list(std::shared_ptr<List<T>> ls) {
                if(nd_ != nullptr)
                        throw std::invalid_argument("YOU CAN T DO THIS\n");
                ls_->concatenate(*ls);
        }
};

namespace Sort {
        template <typename T>
        void merge(T *tab, int b, int m, int e) {
                int s1 = m - b + 1;
                int s2 = e - m;
                T *t1 = new T[s1];
                T *t2 = new T[s2];
                for(int i = 0; i < s1; i++)
                        t1[i] = tab[i + b];
                for(int i = 0; i < s2; i++)
                        t2[i] = tab[i + m + 1];

                int index_t1 = 0;
                int index_t2 = 0;
                int i        = 0;

                for(; index_t1 < s1 && index_t2 < s2; i++) {
                        if(t1[index_t1] < t2[index_t2]) {
                                tab[b + i] = t1[index_t1];
                                ++index_t1;
                                continue;
                        }
                        tab[b + i] = t2[index_t2];
                        ++index_t2;
                }

                for(; index_t1 < s1; index_t1++, i++)
                        tab[b + i] = t1[index_t1];
                for(; index_t2 < s2; index_t2++, i++)
                        tab[b + i] = t2[index_t2];

                delete t1;
                delete t2;
        }

        template <typename T>
        void merge_sort(T *tab, int f, int s) {//first call(tab, 0, tab.size() - 1)
                if(f < s) {
                        int mid = (f + s) >> 1;
                        merge_sort(tab, f, mid);
                        merge_sort(tab, mid + 1, s);
                        merge(tab, f, mid, s);
                }
        }

        template <typename T>
        void insertion_sort(T *a, int s) {
                for(int ii = 1; ii < s; ii++) {
                        int n = 1;
                        int i = ii;
                        while(a[i] < a[ii - n] && n <= ii) {
                                Util::swap_(a[ii - n], a[i]);
                                --i;
                                ++n;
                        }
                }
        }

        template <typename T>
        void insertion_sort(std::vector<T> &a) {//for selection function
                for(int ii = 1; ii < a.size(); ii++) {
                        int n = 1;
                        int i = ii;
                        while(a[i] < a[ii - n] && n <= ii) {
                                Util::swap_(a[ii - n], a[i]);
                                --i;
                                ++n;
                        }
                }
        }

        template <typename T>
        void bubble_sort(T *tab, int n) {
                for(int i = 0; i < n-1; i++) {
                        for(int j = 1; j < n-i; j++) {
                                if(tab[j-1] > tab[j])
                                        Util::swap_(tab[j-1], tab[j]);
                        }
                }
        }

        template <typename T>
        void heap_sort(T *tab, int s) {
                Heap<T> hp = Heap<T>(tab, s);
                for(int i = s-1; i >= 0; i--) {
                        tab[i] = hp.pop();
                }
        }

        template <typename T>
        int partition(T *tab, int p, int r) {
                T pivot = tab[r];
                int n  = 0;
                for(int i = p; i < r; i++) {
                        if(tab[i] > pivot) {
                                ++n;
                        }
                        else{
                                Util::swap_(tab[i - n], tab[i]);
                        }
                }
                Util::swap_(tab[r], tab[r-n]);
                return r-n;
        }

        template <typename T>
        int random_partition(T *tab, int p, int r) {
                srand(getpid());
                auto i = p + (rand() % (p-r));
                Util::swap_(tab[i], tab[r]);
                return partition(tab, p, r);
        }

        template <typename T>
        void quick_sort(T *tab, int b, int e) {//first call(tab, 0, tab.size() - 1)
                if(b < e) {
                        int mid = random_partition(tab, b, e);
                        quick_sort(tab, b, mid-1);
                        quick_sort(tab, mid+1, e);
                }
        }

        template <typename T = uint>
        void counting_sort(T *a, int size_a, int k) {//first call(tab, tab.size(), range of integers)
                uint *tab = new uint[k];
                for(int i = 0; i < size_a; i++) {
                        ++tab[a[i]];
                }
                int index_a = 0;
                for(int i = 0; i < k; i++) {
                        for(int j = 0; j < tab[i]; j++) {
                                a[index_a] = i;
                                ++index_a;
                        }
                }
                delete []tab;
        }
};

bool canPlace(const std::vector<std::vector<int>> &grid, int r, int c)
{
    int i;
    int j;
    auto n = grid.size();

    for(i = 0; i < r; i++)//check the column
        if(grid[i][c] == 1)
            return false;

    for(i = r, j = c; i >= 0 && j < n; i--, j++)//check right-up diagonal
        if(grid[i][j] == 1)
            return false;

    for(i = r, j = c; j >= 0 && i >= 0; i--, j--)//check left-up diagonal
        if(grid[i][j] == 1)
            return false;

    return true;
}

static int solution_n = 0;

bool NQuenns(int n, std::vector<std::vector<int>> &vec, int r) {//vec is a 2D array
        if(r >= n) {
                std::cout << "new : " << ++solution_n << std::endl;
                Util::display(vec);
                return true;
        }

        for(int i = 0; i < n; i++) {
                if(canPlace(vec, r, i)) {
                        vec[r][i] = 1;
                        NQuenns(n, vec, r+1);
                        vec[r][i] = 0;
                }
        }
        return false;
}

template <typename T>
T random_selection(T *a, int b, int e, int i) {
        if(b == e)
                return a[b];

        int med = Sort::random_partition(a, b, e);
        if(med == i)
                return a[i];
        if(i < med)
                return random_selection(a, b, med-1, i);
        return random_selection(a, med+1, e, i);
}

template <typename T>
static std::tuple<int, int, T> max_sub_tab_mid(T *tab, int b, int m, int e) {
        int index_left    = 0;
        auto sum_max_left = std::numeric_limits<T>::min();
        auto sum          = 0;
        for(int i = m; i >= b; i--) {
                sum += tab[i];
                if(sum > sum_max_left) {
                        index_left = i;
                        sum_max_left = sum;
                }
        }
        int index_right   = 0;
        int sum_max_right = std::numeric_limits<T>::min();
        sum               = 0;

        for(int i = m+1; i < e; i++) {
                sum += tab[i];
                if(sum > sum_max_right) {
                        index_right = i;
                        sum_max_right = sum;
                }
        }
        return std::tuple<int, int, T>(index_left, index_right, sum_max_left + sum_max_right);
}

template <typename T>
std::tuple<int, int, T> max_sub_tab(T *tab, int b, int e) {//first call (tab, 0, tab.size())
        if(b == e)
                return std::tuple<int, int, T>(b, e, tab[b]);
        int mid = (b + e) >> 1;
        auto first  = max_sub_tab(tab, b, mid);
        auto second = max_sub_tab(tab, mid + 1, e);
        auto third  = max_sub_tab_mid(tab, b, mid, e);
        if(std::get<1>(first) > std::get<2>(second) && std::get<2>(first) > std::get<2>(third))
                return first;
        if(std::get<2>(second) > std::get<2>(first) && std::get<2>(second) > std::get<2>(third))
                return second;
        return third;
}

int gcd(int a, int b) {
        if(b == 0)
                return a;
        return gcd(b, a % b);
}

std::tuple<int, int, int> extend_euclid(int a, int b) {//first is gcd, second is u third is v
        if(b == 0)
                return std::tuple<int, int, int>(a, 1, 0);
        auto ret = extend_euclid(b, a%b);
        std::swap(std::get<1>(ret), std::get<2>(ret));
        std::get<2>(ret) = std::get<2>(ret) - (a/b)*std::get<1>(ret);
        return ret;
}

int cut_bars(const std::vector<int> &prices, int k) {
        std::vector<long int> tab(k, 0);

        for(int i = 1; i < tab.size(); i++) {
                long int mx = std::numeric_limits<int>::min();
                for(int j = 0; j < i; j++) {
                        mx = Util::max(mx, tab[i-j-1] + prices[j]);
                }
                tab[i] = mx;
        }
        Util::display(tab);
        return tab[k-1];
}
namespace Dynamic {

std::vector<std::pair<std::optional<int>, std::optional<int>>> tab_;//first is without stock and second is with
        int trade(const std::vector<int> &prices, int index, bool have_stock) {
                if(index < 0) {
                        if(have_stock)
                                return std::numeric_limits<short int>::min();
                         return 0;
                }

                int tr1;
                int tr2;
                if(auto opt = tab_[index].first; opt)
                        tr1 = tab_[index].first.value();
                else
                        tab_[index].first = tr1 = trade(prices, index-1, false);

                if(auto opt = tab_[index].second; opt)
                        tr2 = tab_[index].second.value();
                else
                        tab_[index].second = tr2 = trade(prices, index-1, true);

                if(have_stock) {
                        auto buy    = tr1 - prices[index];
                        auto hold   = tr2;
                        return Util::max(buy, hold);
                }
                auto sell   = tr2 + prices[index];
                auto avoid  = tr1;
                return Util::max(sell, avoid);
        }

        int trade(const std::vector<int> &prices) {
                tab_.clear();
                tab_.resize(prices.size());
                return trade(prices, prices.size() - 1, false);
        }

        int trade2(const std::vector<int> &prices) {
                int stock      = std::numeric_limits<short int>::min();
                int cash       = 0;

                for(int i = 0; i < prices.size(); i++) {
                        auto buy    = cash - prices[i];
                        auto hold   = stock;

                        auto sell   = stock + prices[i];
                        auto avoid  = cash;

                        cash  = Util::max(sell, avoid);
                        stock = Util::max(buy, hold);
                }
                return cash;
        }

        std::vector<int> change_making(const std::vector<int> &pieces, int price) {
                std::vector<std::optional<std::vector<int>>> tab(price+1);
                for(int i = 0; i < pieces.size(); i++)
                        tab[pieces[i]] = {pieces[i]};

                for(int i = 0; i < price+1; i++) {
                        auto sz = (tab[i].has_value()) ? tab[i].value().size() : std::numeric_limits<short int>::max();

                        for(int j = 0; j <= (i>>1); j++) {
                                if(!(tab[j].has_value() && tab[i-j].has_value()))
                                        continue;
                                auto tmp = tab[j].value().size() + tab[i-j].value().size();
                                if(tmp > sz)
                                        continue;

                                tab[i] = {std::vector<int>()};
                                tab[i].value().insert(tab[i].value().end(), tab[j].value().begin(), tab[j].value().end());
                                tab[i].value().insert(tab[i].value().end(), tab[i-j].value().begin(), tab[i-j].value().end());
                        }
                }
                return (tab[price].has_value()) ?  tab[price].value() : std::vector<int>();
        }

        int count_change_making_permutations(std::vector<int> &pieces,  int k) {//permutations
                using namespace std;
                if(k <= 0)
                        return 0;
                int ret = 0;
                for(auto e : pieces) {
                        if(e == k) {
                                ++ret;
                                break;
                        }
                        ret += count_change_making_permutations(pieces, k-e);
                }
                return ret;
        }

        int count_change_making_combinations(std::vector<int> &pieces, int k, int index) {//combinasons
                using namespace std;
                if(k == 0)
                        return 1;
                if(index >= pieces.size() || k < 0)
                        return 0;
                int ret = 0;
                for(int i = 0; i <= (k/pieces[index]); i++)
                        ret += count_change_making_combinations(pieces, k-i*pieces[index], index+1);
                return ret;
        }

        int count_expressions(std::vector<int>&v, int t) {
                using namespace std;
                map<int, int> mp;
                mp.insert(std::pair<int, int>(v[0], 1));
                for(int i = 1; i < v.size(); i++) {
                        auto current(mp);
                        mp.clear();
                        for(auto &e : current) {
                                auto sum = e.first + v[i];
                                auto sub = e.first - v[i];
                                if(mp.find(sum) == mp.end())
                                        mp.insert(std::pair<int, int>(sum, e.second));
                                else
                                        mp[sum] += e.second;

                                if(mp.find(sub) == mp.end())
                                        mp.insert(pair<int, int>(sub, e.second));
                                else
                                        mp[sub] += e.second;
                        }
                }

                if(mp.find(t) == mp.end())
                        return 0;
                return mp[t];
        }

        std::vector<std::string> split(const std::string &str, const std::vector<std::string> &dico, std::map<int, std::optional<std::vector<std::string>>> &mp, int index = 0) {
                using namespace std;
                for(int i = 0; i < dico.size(); i++) {
                        int j = 0;
                        auto tmp = index;
                        while(j < dico[i].size() && tmp < str.size() && dico[i][j] == str[tmp]) {
                                ++j;
                                ++tmp;
                        }
                        if(j != dico[i].size())
                                continue;

                        if(tmp == str.size()) {
                                vector<string>r(1, dico[i]);
                                mp.insert(pair<int, optional<vector<string>>>(index, {r}));
                                return r;
                        }

                        if(mp.find(tmp) == mp.end()) {
                                auto vec = split(str, dico, mp, tmp);
                                if(vec.empty())
                                        mp.insert(pair<int, optional<vector<string>>>(tmp, nullopt));
                                else
                                        mp.insert(pair<int, optional<vector<string>>>(tmp, {vec}));
                        }
                        if(!mp[tmp].has_value())
                                continue;
                        auto res = mp[tmp].value();
                        res.push_back(dico[i]);
                        return res;
                }
                return vector<string>();
        }

        int count_BST(int n, std::vector<std::pair<int, int>> &mp) {
                if(n == 0)
                        return 1;
                int ret = 0;
                for(int i = 0; i < n; i++) {
                        if(mp[n-i-1].second < 0)
                                mp[n-i-1].first = count_BST(n-i-1, mp);
                        if(mp[i].second < 0)
                                mp[i].first     = count_BST(i, mp);
                        ret += mp[n-i-1].first * mp[i].first;
                }
                return ret;
        }

        int max_subarray(std::vector<int> &vec) {
                int gain    = 0;
                int penalty = 0;
                int beg     = 0;
                int end     = 0;
                int max     = 0;
                int past    = 0;
                int i       = 0;

                while(true) {
                        if(i == vec.size())
                                break;
                        if(past < 0)
                                gain = 0;
                        else
                                gain = past;

                        while(i < vec.size() && vec[i] >= 0) {
                                gain += vec[i];
                                ++i;
                        }
                        if(gain > max)
                                max = gain;
                        penalty = 0;
                        while(i < vec.size() && vec[i] < 0) {
                                penalty -= vec[i];
                                ++i;
                        }
                        past = gain - penalty;
                }
                return max;
        }

        int max_subarray2(std::vector<int> &vec) {//we assume that therer is no negative numbers
                int max = 1;
                int tmp_min = 1;
                int tmp_max = 1;
                for(int i = 0; i < vec.size() ;i++) {
                        std::cout << "current max = " << max << " and element = " << vec[i] << std::endl;
                        if(vec[i] < 0) {
                                auto l = vec[i] * tmp_min;
                                if(l > max)
                                        max = l;
                                auto k = vec[i] * tmp_max;
                                tmp_min = (k < vec[i]) ? k : vec[i];
                                tmp_max = l;
                        }else if(vec[i] == 0) {
                                tmp_min = 1;
                                tmp_max = 1;
                        }else{
                                tmp_min *= vec[i];
                                tmp_max *= vec[i];
                                if(tmp_max > max)
                                        max = tmp_max;
                        }
                }
                return max;
        }

        std::pair<int, bool> size_palindrom(const std::string & str, int mid, int n) {//[b, e]b and e are include
                int k_o   = 1;
                int k_e   = 0;
                bool even = true;
                bool odd  = true;
                for(int j = 1; j < n && j < mid; j++) {
                        if(odd) {
                                if(str[mid+j] != str[mid-j]) {
                                        odd = false;
                                }else
                                        k_o += 2;
                        }
                        if(even) {
                                if(str[mid-j+1] != str[mid+j])
                                        even = false;
                                else
                                        k_e +=2;
                        }
                        if(!even && !odd)
                                break;
                }
                if(k_e > k_o)
                        return std::pair<int, bool>(k_e, true);
                return std::pair<int, bool>(k_o, false);
        }

        std::string palindromic(const std::string &str) {//in O(n**2)
                using namespace std;
                int m = 0;
                int n = str.size();
                string ret;
                for(int i = 1; i < n; i++) {
                        auto s = size_palindrom(str, i, n);
                        if(s.first > m) {
                                m = s.first;
                                if(s.second) {//it's even
                                        char *tmp =(char*) alloca(m);
                                        for(int k = 0; k < m; k++)
                                                tmp[k] = str[k + i-(s.first>>1)+1];
                                        ret.assign(tmp, m);
                                }else {
                                        char *tmp = (char*)alloca(m);
                                        for(int k = 0; k < m; k++)
                                                tmp[k] = str[k + i-(s.first>>1)];
                                        ret.assign(tmp, m);
                                }
                        }
                }
                return ret;
        }

        std::string manacher(const std::string &str) {
                using namespace std;
                int n     = str.size();
                int r     = 0;
                int e     = 0;
                int mmax  = 0;
                int index = 0;

                for(int i = 1; i < n; i++) {
                        int tmp;
                        if(i > e+r)
                                tmp  = 0;
                        else
                                tmp  = e+r-i;

                        while(tmp+i < n && i-tmp >= 0 && str[tmp+i] == str[i-tmp])
                                ++tmp;

                        auto l =2*tmp-1;
                        if(mmax < l) {
                                index = i-tmp+1;
                                mmax  = l;
                        }
                        if(i+tmp > e+r) {
                                e = i;
                                r = tmp;
                        }
                }
                r    = 0;
                e    = 0;
                for(int i = 1; i < n; i++) {
                        int tmp;
                        if(i > e+r)
                                tmp  = 0;
                        else
                                tmp  = e+r-i;
                        int tmp2 = tmp+1;

                        while(i-tmp >= 0 && i+tmp2 < n && str[i-tmp] == str[i+tmp2]) {
                                ++tmp;
                                ++tmp2;
                        }
                        auto l =2*tmp;
                        if(l > mmax) {
                                mmax  = l;
                                index = i-tmp+1;
                        }
                        if(i+tmp > e+r) {
                                e = i;
                                r = tmp;
                        }
                }
                auto s  = (char*)alloca(mmax);
                for(int i = index; i < index + mmax; ++i)
                        s[i-index] = str[i];
                return string(s, mmax);
        }

        std::string longest_parenthesis(std::string str) {//{{}
                using namespace std;
                vector<int> vec;

                while(str.size() > 0 && str.back() == '{')
                        str.pop_back();

                for(int i = str.size() - 1; i > 0; i--) {
                        if(str[i] == '}' && str[i-1] == '{')
                                vec.push_back(i);
                }

                vector<int> rest(vec.size(), 0);
                for(int i = 0; i < vec.size(); i++) {//we fill in the vec.size()-i-1 case
                        int tmp = vec[i]+1;
                        int cpt = 0;
                        int m   = (i == 0) ? str.size() : vec[i-1]-1;
                        while(tmp < m) {
                                if(str[tmp] == '{')
                                        ++cpt;
                                else if(str[tmp] == '}')
                                        --cpt;
                                ++tmp;
                        }
                        cout << "TOUR : " << i << " AND CPT = " << cpt << endl;
                        if(i == 0)
                                rest.back() = cpt;
                        else
                                rest[vec.size()-i-1] = rest[vec.size()-i] + cpt;
                }

                return string();

        }

};

namespace Graph {
        std::vector<std::list<int>> matrix_to_list(const std::vector<std::vector<int>> &gr) {//wa assume that the graph is oriented
                std::vector<std::list<int>> ret(gr.size());
                for(int i = 0; i < gr.size(); i++)
                        for(int j = 0; j < gr.size(); j++)
                                if(gr[i][j] != 0)
                                        ret[i].push_back(j+1);

                return ret;
        }

        void to_DOT(const std::vector<std::list<int>> &gr, const std::string &name = "graphs.gv") {
                std::ofstream fl(name);
                fl << "digraph G{" <<std::endl;
                for(int i = 0; i < gr.size(); i++) {
                        if(gr[i].empty())
                                continue;

                        fl << "\t"<<i+1<<" -> {";

                        for(auto it = gr[i].begin(); it != gr[i].end(); it++) {
                                if(std::next(it) == gr[i].end())
                                        break;
                                fl << *it << ", ";
                        }
                        fl << gr[i].back() << "};"<<std::endl;
                }
                fl << '}'<<std::endl;
        }

        void to_DOT(const std::vector<std::vector<int>> &gr, const std::string &name = "graphs.gv") {
                auto tmp = matrix_to_list(gr);
                to_DOT(tmp, name);
        }

        std::vector<List<std::pair<int, int>>> to_list(const std::vector<std::list<std::pair<int, int>>> &gr) {
                std::vector<List<std::pair<int, int>>> ret(gr.size());
                for(int i = 0; i < gr.size(); i++) {
                        for(auto &e : gr[i]) {
                                ret[i].push_back(e);
                        }
                }
                return ret;
        }

        std::vector<std::list<std::pair<int, int>>> matrix_to_list_weight(const std::vector<std::vector<int>> &gr) {//wa assume that the graph is oriented
                std::vector<std::list<std::pair<int, int>>> ret(gr.size());

                for(int i = 0; i < gr.size(); i++) {
                        for(int j = 0; j < gr.size(); j++) {
                                if(gr[i][j] != 0)
                                        ret[i].push_back(std::pair<int, int>(j+1, gr[i][j]));
                        }
                }
                return ret;
        }

        std::vector<std::list<std::pair<int, int>>> to_weight(const std::vector<std::list<int>> &gr) {
                std::vector<std::list<std::pair<int, int>>> ret(gr.size());
                for(int i = 0; i < gr.size(); i++) {
                        for(auto e : gr[i]) {
                                ret[i].push_back(std::pair<int, int>(e, 1));
                        }
                }
                return ret;
        }

        std::vector<std::vector<int>> list_to_matrix(const std::vector<std::list<int>> &gr) {
                std::vector<std::vector<int>> ret(gr.size(), std::vector<int>(gr.size()));

                for(int i = 0; i < gr.size(); i++) {
                        for(auto e : gr[i]) {
                                ret[i][e-1] = 1;
                        }
                }
                return ret;
        }

        int outgoing_degree(const std::vector<std::vector<int>> &gr, int s) {
                --s;
                int ret = 0;
                for(int i = 0; i < gr.size(); i++) {
                        if(gr[s][i])
                                ++ret;
                }
                return ret;
        }

        int outgoing_degree(const std::vector<std::list<int>> &gr, int s) {
                int ret = 0;
                for(auto e : gr[s-1]) {
                        ++ret;
                }
                return ret;
        }

        int incoming_degree(const std::vector<std::vector<int>> &gr, int s) {
                --s;
                int ret = 0;
                for(int i = 0; i < gr.size(); i++) {
                        if(gr[i][s] == 1)
                                ++ret;
                }
                return ret;
        }

        int incoming_degree(const std::vector<std::list<int>> &gr, int s) {
                int ret = 0;
                for(int i = 0; i < gr.size(); i++) {
                        for(auto e : gr[i])
                                if(e == s)
                                        ++ret;
                }
                return ret;
        }

        std::vector<std::vector<int>> reverse(std::vector<std::vector<int>> gr) {
                for(int i = 0; i < gr.size(); i++) {
                        for(int j = 0; j < i; j++) {
                                auto tmp = gr[i][j];
                                gr[i][j] = gr[j][i];
                                gr[j][i] = tmp;
                        }
                }
                return gr;
        }

        std::vector<std::list<int>> reverse(const std::vector<std::list<int>> &gr) {
                std::vector<std::list<int>>ret(gr.size());

                for(int i = 0; i < gr.size(); i++) {
                        for(int e : gr[i]) {
                                ret[e-1].push_back(i+1);
                        }
                }
                return ret;
        }

        int is_star(std::vector<std::vector<int>> &vec) {
                std::vector<bool>res(vec.size(), true);
                int j = 0;

                for(int i = 0; i < vec.size(); i++) {
                        if(i == j)
                                continue;
                        if(vec[j][i] == 1) {
                                res[j] = false;
                                ++j;
                        }else
                                res[i] = false;
                }
                for(int i = 0; i < vec.size(); i++)
                        if(res[i])
                                return i+1;
                return 0;
        }

        std::vector<int> BFS(const std::vector<std::list<int>> &gr, int s) {//white = 0, grey = 1, black = 2
                if(s > gr.size())
                        return std::vector<int>();
                std::vector<std::pair<short, short>> attr(gr.size(), std::pair<short, short>{0, std::numeric_limits<short>::min()});//first is color and second is previous
                std::vector<int> ret(gr.size(), std::numeric_limits<short>::max());
                ret[s-1] = 0;

                std::queue<int> fifo;
                fifo.push(s-1);

                while(!fifo.empty()) {
                        int index = fifo.front();
                        fifo.pop();
                        for(auto e : gr[index]) {
                                if(attr[e-1].first == 0) {
                                        attr[e-1].second = index;
                                        attr[e-1].first  = 1;
                                        ret [e-1]        = ret[index] + 1;
                                        fifo.push(e-1);
                                }
                        }
                        attr[index].first = 2;
                }
                return ret;
        }

        void print_path(const std::vector<std::list<int>> &gr, int s, int f) {
                if(s > gr.size() || f > gr.size())
                        return;
                std::vector<std::pair<short, short>> attr(gr.size(), std::pair<short, short>{0, std::numeric_limits<short>::min()});//first is color and second is previous
                std::vector<int> ret(gr.size(), std::numeric_limits<short>::max());
                ret[s-1] = 0;

                std::queue<int> fifo;
                fifo.push(s-1);

                while(!fifo.empty()) {
                        int index = fifo.front();
                        fifo.pop();
                        for(auto e : gr[index]) {
                                if(attr[e-1].first == 0) {
                                        attr[e-1].second = index;
                                        attr[e-1].first  = 1;
                                        ret [e-1]        = ret[index] + 1;
                                        fifo.push(e-1);
                                }
                        }
                        attr[index].first = 2;
                }

                --f;
                --s;

                while(f != s) {
                        std::cout << f + 1<< " ";
                        if(attr[f].second == std::numeric_limits<short>::min())
                                std::cout << "NO PATH EXIST" << std::endl;
                        f = attr[f].second;
                }
                std::cout << std::endl;
        }
static int date_;
        void rec_dfs(const std::vector<std::list<int>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s) {
                for(auto e : gr[s]) {
                        if(color[e-1] == 0) {
                                date[e-1].first  = ++date_;
                                color[e-1]       = 1;
                                prev[e-1]        = s-1;
                                rec_dfs(gr, color, date, prev, e-1);
                        }
                }
                color[s]       = 2;
                date[s].second = ++date_;
        }

        void dfs(const std::vector<std::list<int>> &gr) {
                date_ = 0;
                std::vector<short> color(gr.size());
                std::vector<std::pair<short, short>> date(gr.size());
                std::vector<short>prev(gr.size());


                for(int i = 0; i < gr.size(); i++) {

                        if(color[i] == 0) {
                                date[i].first = ++date_;
                                color[i] = 1;
                                rec_dfs(gr, color, date, prev, i);
                        }
                }
                Util::display(date);
        }

        void rec_dfs_connex(const std::vector<std::list<int>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s) {
                for(auto e : gr[s]) {
                        if(color[e-1] == 0) {
                                date[e-1].first  = ++date_;
                                color[e-1]       = 1;
                                prev[e-1]        = s-1;
                                rec_dfs_connex(gr, color, date, prev, e-1);
                                std::cout <<" "<< e << " ";
                        }
                }
                color[s]       = 2;
                date[s].second = ++date_;
        }

        void connex(const std::vector<std::list<int>> &gr) {
                date_ = 0;
                std::vector<short> color(gr.size());
                std::vector<std::pair<short, short>> date(gr.size());
                std::vector<short>prev(gr.size());


                for(int i = 0; i < gr.size(); i++) {
                        if(color[i] == 0) {
                                date[i].first = ++date_;
                                color[i] = 1;
                                rec_dfs(gr, color, date, prev, i);
                        }
                }

                date_ = 0;
                std::vector<short> color2(gr.size());
                std::vector<std::pair<short, short>> date2(gr.size());
                std::vector<short>prev2(gr.size());

                auto grt = reverse(gr);
                std::vector<int>swap_vec(date.size());
                for(int i = 0; i < date.size(); i++) {
                        auto it    = std::max_element(date.begin(), date.end(), [](const std::pair<short, short> &a, const std::pair<short, short> &b) {return a.second < b.second;});
                        auto index = std::distance(date.begin(), it);
                        it->second = 0;
                        swap_vec[i] = index;
                }

                for(int i = 0; i < grt.size(); i++) {
                        if(color2[swap_vec[i]] == 0) {
                                date2[swap_vec[i]].first = ++date_;
                                color2[swap_vec[i]] = 1;
                                std::cout << "NEW VERTICE : " << swap_vec[i]+1<< " ";
                                rec_dfs_connex(grt,  color2, date2, prev2, swap_vec[i]);
                                std::cout << std::endl;
                        }
                }
        }

        void rec_dfs_topologic_sort(const std::vector<std::list<int>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s, std::vector<int>&tp) {
                for(auto e : gr[s]) {
                        if(color[e-1] == 0) {
                                date[e-1].first  = ++date_;
                                color[e-1]       = 1;
                                prev[e-1]        = s-1;
                                rec_dfs_topologic_sort(gr, color, date, prev, e-1, tp);
                        }
                }
                tp.push_back(s+1);
                color[s]       = 2;
                date[s].second = ++date_;
        }

        std::vector<int> topologic_sort(const std::vector<std::list<int>> &gr) {
                date_ = 0;
                std::vector<short> color(gr.size());
                std::vector<std::pair<short, short>> date(gr.size());
                std::vector<short>prev(gr.size());
                std::vector<int>ret;

                for(int i = 0; i < gr.size(); i++) {
                        if(color[i] == 0) {
                                date[i].first = ++date_;
                                color[i]      = 1;
                                rec_dfs_topologic_sort(gr, color, date, prev, i, ret);
                        }
                }

                std::reverse(ret.begin(), ret.end());
                return ret;
        }

        void rec_dfs_topologic_sort(const std::vector<std::list<std::pair<int, int>>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s, std::vector<int>&tp) {
                for(auto e : gr[s]) {
                        if(color[e.first-1] == 0) {
                                date[e.first-1].first  = ++date_;
                                color[e.first-1]       = 1;
                                prev[e.first-1]        = s-1;
                                rec_dfs_topologic_sort(gr, color, date, prev, e.first-1, tp);
                        }
                }
                tp.push_back(s+1);
                color[s]       = 2;
                date[s].second = ++date_;
        }

        std::vector<int> topologic_sort(const std::vector<std::list<std::pair<int, int>>> &gr) {
                date_ = 0;
                std::vector<short> color(gr.size());
                std::vector<std::pair<short, short>> date(gr.size());
                std::vector<short>prev(gr.size());
                std::vector<int>ret;


                for(int i = 0; i < gr.size(); i++) {
                        if(color[i] == 0) {
                                date[i].first = ++date_;
                                color[i] = 1;
                                rec_dfs_topologic_sort(gr, color, date, prev, i, ret);
                        }
                }

                std::reverse(ret.begin(), ret.end());
                return ret;
        }

        int num_paths(std::vector<std::list<int>> &vec, int s1, int s2) {//only on dag
                --s1;
                int ret = 0;
                for(auto e : vec[s1]) {
                        if(e == s2)
                                ++ret;
                        else
                                ret += num_paths(vec, e, s2);
                }
                return ret;
        }

        std::vector<std::list<std::pair<int, int>>> kruskal(const std::vector<std::vector<int>> &gr) {
                std::vector<std::list<std::pair<int, int>>> ret(gr.size());
                std::vector<std::set<int>> ensemble(gr.size());
                for(int i = 0; i < gr.size(); i++) {
                        ensemble[i].insert(i);
                }

                std::vector<Util::triplet>swap_vec;
                for(int i = 0; i < gr.size(); i++) {
                        for(int j = 0; j < i; j++) {
                                if(gr[i][j] != 0) {
                                        swap_vec.push_back({gr[i][j], i, j});
                                }
                        }
                }
                std::sort(swap_vec.begin(), swap_vec.end(), []( Util::triplet a,  Util::triplet b) {return a.x < b.x;});
                for(auto e : swap_vec) {
                        std::cout << e.x << " " << e.y << " " << e.z << std::endl;
                }

                for(auto e : swap_vec) {
                        if(!Util::is_in(ensemble, e.y, e.z)) {
                                ret[e.y].push_back(std::pair<int, int>(e.x, e.z));
                                ret[e.z].push_back(std::pair<int, int>(e.x, e.y));
                                ensemble[e.z].insert(ensemble[e.y].begin(), ensemble[e.y].end());

                                for(auto it = ensemble[e.z].begin(); it != ensemble[e.z].end(); it++)
                                        ensemble[*it] = ensemble[e.z];
                        }

                }
                return ret;
        }

        std::vector<long int> bellman_ford(const std::vector<std::list<std::pair<int, int>>> &gr, int p) {
                /*
                std::vector<std::vector<int>> graph{{0,6,7,0,0},
                                                    {0,0,8,5,-4},
                                                    {0,0,0,-3,9},
                                                    {0,-2,0,0,0},
                                                    {2,0,0,7,0}};
                */
                --p;
                std::vector<long int> ret(gr.size(), std::numeric_limits<int>::max());
                std::vector<int> prec(gr.size(), std::numeric_limits<int>::max());
                ret[p] = 0;
                for(int j = 0; j < gr.size(); j++) {
                        for(int i = 0; i < gr.size(); i++) {
                                for(auto a : gr[i]) {
                                        ret[a.first - 1] = Util::min(ret[a.first - 1], ret[i] + a.second);
                                }
                        }
                }
                return ret;
        }

        std::vector<long int> bellman_ford_dag(const std::vector<std::list<std::pair<int, int>>> &gr, int p) {
                --p;
                std::vector<long int> ret(gr.size(), std::numeric_limits<int>::max());
                std::vector<int> prec(gr.size(), std::numeric_limits<int>::max());
                auto order = topologic_sort(gr);
                ret[p] = 0;
                for(int i = 0; i < order.size(); i++) {
                        for(auto a : gr[order[i]-1]) {
                                ret[a.first - 1] = Util::min(ret[a.first - 1], ret[order[i]-1] + a.second);
                        }

                }
                return ret;
        }

        std::vector<long int> dikjstra(const std::vector<std::list<std::pair<int, int>>> &gr, int p) {
                --p;
                std::vector<long int> ret(gr.size(), std::numeric_limits<int>::max());
                Fibonacci_heap<std::pair<int, long int>> hp(Util::pair_comparaison_dikjstra<int, long int>);
                std::vector<Fibonacci_root< std::pair<int, long int> >*> fd(gr.size());
                std::vector<bool> erase(gr.size(), false);
                ret[p]        = 0;

                for(int i = 0; i < gr.size(); i++) {
                        fd[i] = hp.insert(std::pair<int, long int>(i, ret[i]));
                }
                bool ok = true;

                while(!hp.empty()) {
                        auto s         = hp.get_min()->get_value();
                        fd[s.first]    = nullptr;
                        erase[s.first] = true;
                        hp.pop();

                        for(auto &e : gr[s.first]) {
                                if(ret[e.first-1] > ret[s.first] + e.second) {
                                        ret[e.first-1] = ret[s.first] + e.second;
                                        if(!erase[e.first-1])
                                                hp.reduce(fd[e.first - 1], std::pair<int, long int>(e.first-1, ret[s.first] + e.second));
                                }
                                ok = false;
                        }
                }
                return ret;
        }

        std::vector<std::vector<long int>> pow_paths(std::vector<std::vector<long int>> gr1, const std::vector<std::vector<long int>> &gr2) {
                for(int l = 0; l < gr1.size(); l++) {
                        for(int i = 0; i < gr1.size(); i++) {
                                long int tmp = std::numeric_limits<int>::max();
                                for(int j = 0; j < gr1.size(); j++)
                                        tmp = std::min(tmp, gr1[l][j]+gr2[j][i]);
                                gr1[l][i] = tmp;
                        }
                }
                return gr1;
        }

        std::vector<std::vector<long int>> shortest_paths_pairs_naive(std::vector<std::vector<long int>> gr) {
                for(int i = 0; i < gr.size(); i++)
                        for(int j = 0; j < gr.size(); j++)
                                if(i != j && gr[i][j] == 0)
                                        gr[i][j] = std::numeric_limits<int>::max();

                auto n  = gr.size();
                auto bs = gr;
                while(n > 0) {
                        if(n % 2 == 1)
                                gr = pow_paths(gr, bs);

                        gr = pow_paths(gr, gr);
                        n>>=1;
                }
                return gr;
        }

        std::vector<std::vector<long int>> floyd_warshall(const std::vector<std::vector<long int>> &gr) {
                auto past(gr);
                for(int i = 0; i < gr.size(); i++)
                        for(int j = 0; j < gr.size(); j++)
                                if(i != j && past[i][j] == 0)
                                        past[i][j] = std::numeric_limits<int>::max();

                std::vector<std::vector<long int>> pred(gr);
                for(int i = 0; i < pred.size(); i++) {
                        for(int j= 0; j < pred.size(); j++) {
                                if(pred[i][j] != 0)
                                        pred[i][j] = i+1;
                        }
                }
                auto n(past);
                for(int k = 0; k < gr.size(); k++) {
                        n = past;
                        for(int i = 0; i < n.size(); i++) {
                                for(int j = 0; j < n.size(); j++) {
                                        if(n[i][j] > past[i][k] + past[k][j]) {
                                                pred[i][j] = pred[k][j];//compute predecessor matrice
                                                n[i][j]    = past[i][k] + past[k][j];
                                        }
                                }
                        }
                        past = n;
                }
                //Util::display(pred);
                //std::cout << std::endl;

                return n;
        }

        std::vector<std::vector<int>> transitive_closure(const std::vector<std::vector<int>> &gr) {
                auto past(gr);
                auto n(past);

                for(int i = 0; i < gr.size(); i++)
                        past[i][i] = 1;

                for(int k = 0; k < gr.size(); k++) {
                        n = past;
                        for(int i = 0; i < n.size(); i++) {
                                for(int j = 0; j < n.size(); j++) {
                                        n[i][j] = past[i][j] || (past[i][k] && past[k][j]);
                                }
                        }
                        past = n;
                }
                return n;
        }

        std::vector<std::list<std::pair<int, int>>> johnson(std::vector<std::list<std::pair<int, int>>> gr) {
                std::list<std::pair<int, int>> pb;
                for(int i = 0; i < gr.size(); i++)
                        pb.push_back(std::pair<int, int>(i+1, 0));
                gr.push_back(pb);

                auto res = bellman_ford(gr, gr.size());

                for(int i = 0; i < gr.size() - 1; i++)
                        for(auto &e : gr[i])
                                e.second = e.second + res[i] - res[e.first-1];

                std::vector<std::list<std::pair<int, int>>> ret(gr.size()-1);
                for(int i = 0; i < gr.size()-1; i++) {
                        auto it = dikjstra(gr, i+1);
                        for(int j = 0; j < it.size()-1; j++)
                                ret[i].push_back(std::pair<int, int>(j, it[j]+res[j]-res[i]));
                }
                return ret;
        }

        bool rec_dfs_DAG(const std::vector<std::list<std::pair<int, int>>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s) {
                auto ret = false;
                for(auto e : gr[s]) {
                        if(color[e.first - 1] == 1) {
                                std::cout << "FROM HERE : " << s +1 << "AND : " << e.first << std::endl;
                                ret = true;
                        }
                        if(color[e.first-1] == 0) {
                                date[e.first-1].first  = ++date_;
                                color[e.first-1]       = 1;
                                prev[e.first-1]        = s-1;
                                ret |= rec_dfs_DAG(gr, color, date, prev, e.first-1);
                        }
                }
                color[s]       = 2;
                date[s].second = ++date_;
                return ret;
        }


        bool is_DAG(const std::vector<std::vector<int>> &gr) {
                auto ls = matrix_to_list_weight(gr);

                std::vector<short> color(gr.size());
                std::vector<std::pair<short, short>> date(gr.size());
                std::vector<short>prev(gr.size());
                bool ret = false;

                for(int i = 0; i < gr.size(); i++) {
                        if(color[i] == 0) {
                                date[i].first = ++date_;
                                color[i] = 1;
                                ret |= rec_dfs_DAG(ls, color, date, prev, i);
                        }
                }
                return !ret;

        }

        void rec_dfs_transitive_closure(const std::vector<std::list<std::pair<int, int>>> &gr, std::vector<short> &color, std::vector<std::pair<short, short>> &date, std::vector<short> &prev, int s, std::vector<std::list<int>>&tp) {
                for(auto &e : gr[s]) {
                        if(color[e.first-1] == 0) {
                                date[e.first-1].first  = ++date_;
                                color[e.first-1]       = 1;
                                prev[e.first-1]        = s-1;
                                rec_dfs_transitive_closure(gr, color, date, prev, e.first-1, tp);
                        }
                        tp[s].insert(tp[s].end(), tp[e.first-1].begin(), tp[e.first-1].end());
                }
                color[s]       = 2;
                date[s].second = ++date_;
        }

        std::vector<std::list<int>> random_dag(int max_levels = 5, int max_per_level = 5) {
                srand(getpid());
                auto initial = max_per_level;
                int past_size = 5;
                std::vector<std::list<int>>ret(past_size);

                for(int i = 0; i < max_levels; i++) {
                        auto actual_size = ret.size();
                        ret.resize(ret.size() + max_per_level);

                        for(int j = 0; j < past_size; j++) {//on itere sur les elements du rang d'avant
                                auto tour = 3;

                                for(int k = 0; k < tour; k++) {//on itere donc sur tour

                                        auto rd = rand() % max_per_level;//on choisit un element aleatoirement dans le rang suivant si il n'est pas deja dans la liste
                                        if(!Util::is_in<int>(ret[actual_size - j - 1], rd + actual_size))
                                                ret[actual_size - j-1].push_back(rd + actual_size);
                                }
                        }

                }
                return ret;
        }

        std::vector<std::list<std::pair<int, int>>> ford_fulkerson(const std::vector<std::list<std::pair<int, int>>> &gr) {//the flot must travel from the first to the last element
                std::vector<std::list<std::pair<int, int>>> residual(gr);
                std::vector<std::list<std::pair<int, int>>> ret(gr.size());
                for(int i = 0; i < gr.size(); i++) {
                        for(auto e : gr[i]) {
                                residual[e.first-1].push_back(std::pair<int, int>(i+1, 0));
                                ret[i].push_back(std::pair<int, int>(e.first, 0));
                        }
                }

                std::queue<int>fifo;
                std::vector<int> color(gr.size());
                std::vector<int> pred(gr.size(), 0);

                while(true) {
                        for(int i = 0; i < gr.size(); i++) {
                                pred[i] = -1;
                                color[i] = 0;
                        }

                        fifo.push(0);
                        color[0] = 2;
                        while(!fifo.empty()) {
                                auto e = fifo.front();
                                fifo.pop();
                                for(auto &ele : residual[e]) {
                                        if(color[ele.first-1] == 0 && ele.second > 0) {
                                                fifo.push(ele.first-1);
                                                color[ele.first-1] = 2;
                                                pred[ele.first-1]  = e;
                                        }
                                }
                        }
                        int n  = pred.back();
                        int af = pred.size();
                        if(n == -1)
                                break;

                        int residual_capacity = std::numeric_limits<int>::max();
                        while(n != -1) {
                                for(auto e : residual[n]) {
                                        if(e.first == af) {
                                                residual_capacity = (residual_capacity > e.second)  ? e.second : residual_capacity;
                                                af = n+1;
                                                n  = pred[n];
                                                break;
                                        }

                                }
                        }

                        n  = pred.back();
                        af = pred.size();
                        while(n != -1)  {//flot is improved
                                bool f = false;
                                for(auto &e : ret[n]) {
                                        if(e.first == af) {
                                                e.second += residual_capacity;
                                                af = n+1;
                                                n  = pred[n];
                                                f = true;
                                                break;
                                        }
                                }
                                if(f)
                                        continue;
                                for(auto &e : ret[af-1]) {
                                        if(e.first == n + 1) {
                                                e.second -= residual_capacity;
                                                af = n+1;
                                                n  = pred[n];
                                                break;
                                        }
                                }
                        }
                        //now we have to update the residual flot now

                        n  = pred.back();
                        af = pred.size();
                        while(n != -1)  {//flot is improved
                                for(auto &e : residual[n]) {
                                        if(e.first == af) {
                                                e.second -= residual_capacity;
                                                break;
                                        }
                                }
                                for(auto &e : residual[af-1]) {
                                        if(e.first == n + 1) {
                                                e.second += residual_capacity;
                                                break;
                                        }
                                }
                                af = n+1;
                                n  = pred[n];
                        }
                }
                return ret;
        }

};

#define N 40

int main() {
        std::string str("{{}");
        std::cout << "MAX PALINDROM = " << Dynamic::longest_parenthesis(str) << std::endl;
        return 0;
}
