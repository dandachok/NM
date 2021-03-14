#pragma once

#include "MatrixAlgo.hpp"

class TPoint {
    TPoint (double x, double y) :x(x), y(y) {}

private:
    double x;
    double y;
};

using vpoint = std::vector<TPoint>;


vdouble GetY(vdouble x, double (*f)(double)) {
    vdouble y(x.size());
    for (int i = 0; i < x.size(); ++i) {
        y[i] = f(x[i]);
    }
    return y;
}

class Diff {
public:
    Diff(vdouble x, vdouble y): x(x), y(y) {}

    double dif(double x0) {
        int n = x.size();
        int i = getInterval(x0);
        double d = dy(i) / dx(i);
        d += (dydx(i + 1) - dydx(i)) * (2*x0 - x[i] - x[i+1]) / (x[i + 2] - x[i]);
        return d;
    }

    double dif2(double x0) {
        int n = x.size();
        int i = getInterval(x0);
        double d;
        d = 2*(dydx(i + 1) - dydx(i)) / (x[i + 2] - x[i]);
        return d;
    }

private:

    double dx(int i) {
        return x[i + 1] - x[i];
    }
    
    double dy(int i) {
        return y[i + 1] - y[i];
    }
        
    double dydx(int i) {
        return dy(i) / dx(i);
    }

    int getInterval(double x0) {
        int n = x.size();
        int i;
        for (i = 0; i < n - 1; ++i) {
            if (x[i] <= x0 && x0 <= x[i + 1]) {
                break;
            }
        }
        return i;
    }
    
    vdouble x;
    vdouble y;
};

class Integrate {
public:
    Integrate(double (*f)(double)) : f(f) {}
    Integrate(vdouble x, vdouble y) :x(x), y(y){}

    void makeNet(double l, double r, double h) {
        x = GetX(l, r, h);
        y = GetY(x, f);
    }

    double RectangleMethod (double l, double r, double h) {
        makeNet(l, r, h);
        int n = x.size();
        double res = 0;
        for (int i = 0; i < n - 1; ++i) {
            double d =  f((x[i] + x[i + 1]) / 2);
            res += d;
        }
        res *= h;
        return res;
    }

    double Trapezium (double l, double r, double h) {
        makeNet(l, r, h);
        int n = x.size();
        double res = Sum(y);
        res -= (y[0] + y[n - 1]) / 2;
        res *= h;
        return res;
    }

    double Simpson (double l, double r, double h) {
        makeNet(l, r, h);
        int n = x.size();
        double res = y[0] + y[n - 1];
        for (int i = 1; i < n - 1; ++i) {
            res += i % 2? 4*y[i]: 2*y[i];
        }
        res *= h / 3;
        return res;
    }

    double RungeRombergOpt (double f1, double f2, double h1, double h2) {
        double k = h2 / h1;
        return f1 + (f1 - f2) / (k*k - 1);
    }
private:

    vdouble GetX(double l, double r, double h) {
        double x0 = l;
        vdouble x;
        while (x0 <= r) {
            x.push_back(x0);
            x0 += h;
        }
        return x;
    }

    double (*f)(double) = nullptr;
    vdouble x;
    vdouble y;
};
