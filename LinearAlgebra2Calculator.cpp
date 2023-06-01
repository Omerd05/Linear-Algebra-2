#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <string>
#include <cstdlib>
#include <ctime>
#include <regex>
#include <istream>
#include <sstream>

using namespace std;
typedef long double Real;
typedef long long Natrual;
#define loop(x, n) for(int x = 0; x < n; ++ x)

const Real epsilon = 1e-6;
const Real delta = 1e-6;
const Real pi = 3.14159265;

struct Polynomial
{
    Natrual deg = 0;
    vector<Real> coef;
    set<Natrual> exist; //This set will maintain and contain which indexs aren't zero. We'll use this set in order to improve our time complexity.
    map<Real, Natrual> roots;
    

    Polynomial() {}

    Polynomial(const Polynomial& p) {
        deg = p.deg;
        coef.resize(p.coef.size(),0);

        for (int i = 0; i < p.coef.size(); i++) {
            coef[i] = p.coef[i];
        }
        for (auto u : p.exist) {
            exist.insert(u);
        }
        for (auto u : p.roots) {
            roots.insert(u);
        }
    }

    //need to fix ASAP
    Polynomial(string poly) { //Excepting to get a function like this - p(x) = 1x^0 + 1.5x^1

        string::iterator end_pos = remove(poly.begin(), poly.end(), ' ');
        poly.erase(end_pos, poly.end());

        vector<Natrual> powers;

        int i = 0;
        while (i++ < poly.length()-1){
            if (poly[i + 1] == '^') {
                powers.push_back(poly[i + 2] - '0');
            }
        }

        deg = 0;
        for (Natrual val : powers)deg = max(deg, val);

        char delimiter = '+';
        std::istringstream iss(poly);
        std::vector<std::string> tokens;

        std::string token;

        Natrual cnt = 0;

        while (getline(iss, token, delimiter)) {
            if(cnt%2!=1)tokens.push_back(token);
        }


        coef.resize(deg + 1, 0);
        for (int i = 0; i < powers.size(); i++) {
            tokens[i].pop_back(); tokens[i].pop_back(); tokens[i].pop_back();
            coef[powers[i]] = stod(tokens[i]);
            exist.insert(powers[i]);
        }
    }

    Real P(Polynomial p, Real x) { //value at point x.
        Real val = 0;
        for (Natrual i : exist) {
            val += coef[i] * pow(x, i);
        }
        return val;
    }

    Polynomial add(Polynomial p, Polynomial q) {    //WORKS ONLY FOR NON-NEGATIVE POLYNOMIALS, NEED TO FIX ASAP
        Polynomial result = Polynomial();
        result.coef.resize(max(p.deg, q.deg) + 1,0);

        set<Natrual> seen;
        for (Natrual i : p.exist){ //(int i = 0; i < coef.size(); i++) {
            result.coef[i] = p.coef[i] + q.coef[i];
            
            if (result.coef[i] != 0) {
                result.exist.insert(i);
                
                result.deg = max(result.deg, i);
            }

            seen.insert(i);
        }

        for (Natrual i : q.exist) {
            if (seen.count(i)) { continue; }

            result.coef[i] = p.coef[i] + q.coef[i];

            if (result.coef[i] != 0) {
                result.exist.insert(i);

                result.deg = max(result.deg, i);
            }

        }

        seen.clear();

        return result;
    }

    Polynomial multi(Polynomial p, Polynomial q) {
        Polynomial result = Polynomial();
        result.deg = p.deg + q.deg;
        result.coef.resize(deg + 1,0);
        
        for (Natrual i : p.exist) {
            for (Natrual j : q.exist) {
                result.coef[i + j] += coef[i] * coef[j];
            }
        }
    }
    
    pair<Polynomial,Polynomial> divison (Polynomial p, Polynomial q) { // Returns {Quotient,Reminder},

        Polynomial quo = Polynomial("0x^0");
        Polynomial rem = Polynomial(p);

        if (q.deg == 0) {
            return { quo,rem };
        }

        while (rem.deg >= q.deg) {
            

            Real t = rem.coef[rem.deg] / q.coef[q.deg];

            auto factor = Polynomial("" + to_string(t) + "x^" + to_string(rem.deg - q.deg));
            
            quo = add(quo,factor);
            
            rem = add(rem, multi(multi(rem, Polynomial("-1x^0")), q));
            
        }
        
        return {quo,rem};
    }

    Polynomial derivative(Polynomial p) {
        Polynomial result = Polynomial();
        result.deg = p.deg - 1;
        result.coef.resize(result.deg+1, 0);
        for (Natrual i : p.exist) {
            if(i>0)result.exist.insert(i - 1);
        }
        for (Natrual i : result.exist) {
            result.coef[i] = p.coef[i + 1] * (i+1);
        }
        for (pair<Real,Natrual> u : p.roots) {
            if (u.second > 0) {
                result.roots.insert({ u.first,u.second - 1 });
            }
        }
        return result;
    }
    
    Real bisection(Polynomial p, Real xLeft, Real xRight) { //Assuming p(xLeft) * p(xRight) < 0.
        Real fxLeft = P(p,xLeft);
        Real fxRight = P(p, xRight);
        Real mid = (xLeft + xRight) / 2;
        Real fmid = P(p, mid);

        Natrual t = 0;

        while (abs(fmid) > delta && t++ < 40) {


            Real mid = (xLeft + xRight) / 2;
            Real fmid = P(p, mid);

            if (abs(fmid) < delta || t == 39) {
                if (abs(mid - ceil(mid)) < epsilon) { mid = ceil(mid); }
                return mid;
            }

            if (fxLeft < 0) {
                if (fmid < 0) {
                    xLeft = mid;
                    fxLeft = fmid;
                }
                else {
                    xRight = mid;
                    fxRight = fmid;
                }
            }
            else {
                if (fmid > 0) {
                    xLeft = mid;
                    fxLeft = fmid;
                }
                else {
                    xRight = mid;
                    fxRight = fmid;
                }
            }

        }
        return 0;
    }

    Real newtonRhapson(Polynomial p) {
        Natrual t = 0;
        Polynomial pd = derivative(p);

        while (t++ < 50) {

            Natrual j = 0;
            Real x = rand() % 5000;
            Real fx = P(p, x);
            Real fdx = P(pd, x);

            while (abs(fx) > delta && j++ < 50) {

                x = x - Real(fx / fdx);
                Real fx = P(p, x);
                Real fdx = P(pd, x);

                if (abs(fx) < delta) {
                    if (abs(x - ceil(x)) < epsilon) { x = ceil(x); }
                    return x;
                }
            }
            
        }
        return 1e9;
    }

    Real rootFinder(Polynomial p) { //we'll assume our polynomial is fixed. If our polynomial is with an odd degree we'll use bisection method otherwise newton-rhapson and pray

        Natrual t = 0;
        if (p.deg % 2) { //When our degree is odd, there is an easy way to find the interval, as when x->inf there exists an x0 st for every x > x0 : f(x)*f(-x) < 0
                         // because when p(x) -x->inf> inf then p(-x) -x->inf> -inf (and the opposite way).
            Real xL = -1;
            Real xR = 1;

            while (P(p, xL) * P(p, xR) > 0 && t++<40) {
                xL *= 2;
                xR *= 2;
            }
            
            return bisection(p, xL, xR);
        }
        else {
            Real x = newtonRhapson(p);
            if (x == 0) {
                if (abs(P(p, x)) >= epsilon) {
                    x = 1e9;
                }
            }
            return x;
        }
    }

    public: void Print() {
        string result = "";
        for (Natrual i : exist) {
            result += to_string(coef[i]) + "x^" + to_string(i);
            if (i < deg) { result += " + "; }
        }
        cout << result << '\n';
    }
};

//Only for square matrix, and above R.

struct Matrix { 
    Natrual n;
    Real det;
    vector<vector<Real>> M; //(n,vector<Real>(n));

    /*Matrix(const Matrix& mat) {
        n = mat.n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = mat.M[i][j];
            }
        }
    }*/

    Matrix(Natrual sz, vector<vector<Real>> tbl) {
        n = sz;
        M.resize(n,vector<Real>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = tbl[i][j];
            }
        }
    }

    void print() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << M[i][j] << " ";
            }
            cout << endl;
        }
    }

    void add(Matrix Q) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] += Q.M[i][j];
            }
        }
    }

    void scalar(Real u) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] *= u;
            }
        }
    }

    Matrix transpose() {
        vector<vector<Real>> T(n,vector<Real>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                T[i][j] = M[j][i];
            }
        }
        return Matrix(n, T);
    }

    Matrix multi(Matrix B) {
        Matrix BT = B.transpose();
        Matrix result = Matrix(n,vector<vector<Real>>(n,vector<Real>(n,0)));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Real sum = 0;
                for (int k = 0; k < n; k++) {
                    sum += M[i][k] + BT.M[j][k];
                }
                result.M[i][j] = sum;
            }
        }
        return result;
    }

    //swapping rows j and i.
    void swapR(int i, int j) { 
        for (int k = 0; k < n; k++) {
            auto temp = M[i][k];
            M[i][k] = M[j][k];
            M[j][k] = temp;
        }
    }

    //R_i <- R_i + u*R_j
    void Rplus(int i, int j, int u) {
        for (int k = 0; k < n; k++) {
            M[i][k] += u*M[j][k];
        }
    }

    Matrix reducedform() {
        Matrix A = Matrix(n, M);
        int i = 0;
        for (int j = 0; j < n; j++){

            int cand = i;
            while (cand < n && A.M[cand][j] == 0)cand++;
            if (cand == n)continue;
            A.swapR(cand, i);

            int runner = i+1;
            while (runner < n) {
                A.Rplus(runner, i, -1 * A.M[runner][j] / A.M[i][j]);
                runner++;
            }

            i++;
        }
        return A;
    }

    void setDet() {
        Matrix A = reducedform();
        Real mul = 1;
        for (int i = 0; i < n; i++) {
            mul *= M[i][i];
        }
        det = mul;
    }
};

int main()
{
    //vector<vector<Real>> tbl = { {1,2,3},{4,5,6},{7,8,8} };
    //Matrix A = Matrix(3, tbl);
    //A.reducedform().print();
    //auto p = Polynomial("-4x^0 + 1x^2");
    //p.Print();
    //cout << p.rootFinder(p) << endl;
}
