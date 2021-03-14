#pragma once

#include <vector>

#include "MatrixAlgo.hpp"

class TPolynom {
public:
    TPolynom(int n);
    TPolynom(vdouble coef);
    TPolynom(vdouble& coef);

    double& operator[] (int i) {
        return coef[i]; 
    }
    const double& operator[] (int i) const {
        return coef[i]; 
    }

    TPolynom operator* (const TPolynom& b) const {
        int n = b.coef.size();
        int m = coef.size();
        TPolynom c(n + m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
            c[i + j] += coef[i] * b[j];
            }
        }
        return c;
    }

private:
    std::vector<double> coef;
};

using polynom = std::vector<double>;


polynom operator* (double v, polynom& b) {
    for (int i = 0; i < b.size(); ++i) {
        b[i] *= v;
    }
    return b;
}

vdouble CoefPolynomLagrange (const vdouble& x, const vdouble& y) {
    int n = x.size();
    vdouble w(n, 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                w[i] *= x[i] - x[j];
            }
        }
        w[i] = y[i] / w[i];
    }
    return w;
}

vdouble CoefPolynomLagrange (vdouble x, double (*f)(double)) {
    int n = x.size();
    vdouble y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = f(x[i]);
    }

    return CoefPolynomLagrange(x, y);
}

double PolynomLagrange(const vdouble& x, const vdouble& y, double x0) {
    int n = x.size();
    vdouble coef = CoefPolynomLagrange(x, y);
    double y0 = 0;
    for (int i = 0; i < n; ++i) {
        double p = 1;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                p *= x0 - x[j];
            }
        }
        y0 += p * coef[i];
    }

    return y0;
}

double PolynomLagrange (const vdouble& x, double (*f)(double), double x0) {
    int n = x.size();
    vdouble y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = f(x[i]);
    }

    return PolynomLagrange(x, y, x0);
}


class PolynomNewton {
public:

    PolynomNewton (const vdouble& x, const vdouble& y) {
        ft = vvdouble(x.size(), vdouble());
        X = x;
        int n = x.size();
        for (int i = 0; i < n; ++i) {
            ft[i].push_back(y[i]);
        }
        for (int j = 0; j < n - 1; ++j) {
            for (int i = 0; i < n - j - 1; ++i) {
                double f = ft[i][j] - ft[i + 1][j];
                f /= x[i] - x[i + j + 1];
                ft[i].push_back(f);
            }
        }
    }

    void AddPoint(double x, double y) {
        int n = X.size();
        X.push_back(x);
        ft.push_back(vdouble());
        ft[n].push_back(y);
        for (int i = n; i > 0; --i) {
            double f = ft[i - 1][n - i] - ft[i][n - i];
            f /= X[i - 1] - X[n];
            ft[i - 1].push_back(f);
        }
    }

    const vvdouble& GetFT() const {
        return ft;
    }

    double operator() (double x0) const {
        int n = X.size();
        double y0 = 0;
        double p = 1;
        for (int i = 0; i < n; ++i) {
            y0 += ft[0][i] * p;
            p *= (x0 - X[i]);

        }
        return y0;
    }
private:
    vvdouble ft;
    vdouble X;
    vdouble Y;
};


class CubicSplines {
public:
    CubicSplines (vdouble x, vdouble y) :
    x(x) {
        int n = x.size();
        c = SweepMethod(GetSweepMatrix(x, y));
        c.insert(c.begin(), 0);
        a = y;
        a.pop_back();
        b.resize(n - 1);
        d.resize(n - 1);
        for (int i = 0; i < n - 2; ++i) {
            double h = (x[i + 1] - x[i]);
            b[i] = (y[i + 1] - y[i]) / h;
            b[i] -=  h * (c[i + 1] + 2*c[i]) / 3;
            d[i] = (c[i + 1] - c[i]) / h / 3;
        }
        double h = x[n - 1] - x[n - 2];
        b[n - 2] = (y[n - 1] - y[n - 2]) / h - 2*h*c[n - 2] / 3;
        d[n - 2] = - c[n - 2] / h / 3;
    }

    double operator() (double x0) {
        double y0;
        int i;
        for (i = 0; i < x.size() - 1; ++i) {
            if (x[i] <= x0 && x0 <= x[i + 1]) {
                break;
            }
        }
        double dx = x0 - x[i];
        y0 = a[i] + b[i]*dx + c[i]*pow(dx, 2) + d[i]*pow(dx, 3);
        return y0;
    }
private:

    vvdouble GetSweepMatrix (const vdouble& x, const vdouble& y) {
        int n = x.size();
        vvdouble m(n - 2, vdouble(4));
        vdouble h = GetH(x);
        m[0][0] = 0;
        m[0][1] = 2 * (h[1] + h[2]);
        m[0][2] = h[2];
        m[0][3] = GetD(y, h, 0);

        for (int i = 1; i < n - 2; ++i) {
            m[i][0] = h[i + 1];
            m[i][1] = 2 * (h[i + 1] + h[i + 2]);
            m[i][2] = h[i + 2];
            m[i][3] = GetD(y, h, i);
        }
        m[n - 3][2] = 0;
        return m;
    }

    vdouble GetH (const vdouble& x) {
        int n = x.size();
        vdouble h(x.size());
        h[0] = 1;
        for (int i = 1; i < n; ++i) {
            h[i] = x[i] - x[i - 1];
        }
        return h;
    }

    double GetD (const vdouble& y, const vdouble& h, int i) {
        double res = (y[i + 2] - y[i + 1]) / h[i + 2];
        res -= (y[i + 1] - y[i]) / h[i + 1];
        return 3*res;
    }

    vdouble x;
    vdouble a;
    vdouble b;
    vdouble c;
    vdouble d;
};


class MNK {
public:
    MNK (vdouble x, vdouble y, int k):
    x(x), y(y), k(k)  {
        vvdouble m = GetNormalMat();
        int n = m.size();
        vdouble b(n);
        for (int i = 0; i < n; ++i) {
            b[i] = m[i][n];
            m[i].pop_back();
        }
        a = GausseMethod(m, b);
        //std::cout << a;
    }

    double F (double x0) {
        int n = a.size();
        double y0 = 0;
        for (int i = 0; i < n; ++i) {
            y0 += a[i]*pow(x0, i);
        }
        return y0;
    }

    double operator() (double x0) {
        return F(x0);
    }

    vdouble operator() (vdouble x) {
        vdouble y;
        for (auto i: x) {
            y.push_back(F(i));
        }
        return y;
    }

    double SSE () {
        int n = x.size();
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += pow(F(x[i]) - y[i], 2);
        }
        return sum;
    }
private:
    vvdouble GetNormalMat() {
        int n = k + 1;
        vdouble coef(2*n - 1);
        vvdouble m(n, vdouble(n));
        for (int i = 0; i < 2*n - 1; ++i) {
            coef[i] = Sum(Pow(x, i));
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                m[i][j] = coef[i + j];
            }
        }

        for (int i = 0; i < n; ++i) {
            m[i].push_back(Sum(Pow(x, i) * y));
        }
        return m;
    }

    vdouble x;
    vdouble y;
    vdouble a;
    int k;  // pow polinom
};
