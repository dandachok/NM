#include <iostream>

#include "NMLibs.hpp"

const double e = 2.71828;

double f (double x, double y, double z) {
    return -(4*x*x + 2)*y - 4*x*z;
}

double f_real(double x) {
    return (1 + x) / std::pow(e, x*x);
}

class MyFunc final: public Function3 {
    public: double operator() (double x, double y, double z) override {
        return (4 * y - 4 * x * z) / (2 * x + 1);
    }
};

class PFun final: public Function1 {
    public: double operator() (double x) {
        return 4 * x / (2 * x + 1);
    }
};

class QFun final: public Function1 {
    public: double operator() (double x) {
        return -4 / (2 * x + 1);
    }
};

class RealF final: public Function1 {
    public: double operator() (double x) override {
        return x + std::pow(e, -2*x);
    }
};

int main() {
    double left = 0;
    double right = 1;
    double h = 0.1;
    double lz0 = -1;
    double a = 1;
    double b = 2;
    double c = 3;
    double eps = 0.1;
    LDE re(a, b, c);

    MyFunc f;
    PFun p;
    QFun q;
    BVProblemL2R3 task(&f, left, right, lz0, re);
    
    ShootingMethodL2R3 s(task);
    s.Calc(h, e);
    std::cout << '\n' << s.GetY() << '\n';

    FDMethod fd(&p, &q, task);
    fd.Calc(h);
    std::cout << '\n' << fd.GetY() << '\n';

    RealF rf;
    vdouble x = GetValInRange(left, right, h);
    vdouble ry = GetFValInRange(x, &rf);
    std::cout << ry;
//
    //EilerMethod2 m(f, left, right, y0, z0);
    //m(h);
    //vdouble y = m.GetY();
    //std::cout << "Eiler method:\n" << y;
//
    //RungeKutta2 rk(f, left, right, y0, z0);
    //rk(h);
    //vdouble rky = rk.GetY();
    //std::cout << "Runge-Kutta method 4 order:\n" << rky;
//
    //AdamsMethod a(f, left, right, y0, z0);
    //a(h);
    //vdouble ay = a.GetY();
    //std::cout << "Adams method:\n" << ay;
//
    //vdouble x = GetValInRange(left, right, h);
    //y = GetY(x, f_real);
    //std::cout << "Real funcrion:\n" << y;



    return 0;
}
