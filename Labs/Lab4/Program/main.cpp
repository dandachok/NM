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
 
 class PFunc final: public Function1 {
     public: double operator() (double x) {
         return 4 * x / (2 * x + 1);
     }
 };
 
 class QFunc final: public Function1 {
     public: double operator() (double x) {
         return -4 / (2 * x + 1);
     }
 };
 
 class RealFunc final: public Function1 {
     public: double operator() (double x) override {
         return x + std::pow(e, -2*x);
     }
 };

//class MyFunc: public Function3 {
//    double operator() (double x, double y, double z) override {
//        return 2*y / (x*x*(x + 1));
//    }
//};
//
//class PFunc: public Function1 {
//    double operator() (double x) override {
//        return 0;
//    }
//};
//
//class QFunc: public Function1 {
//    double operator() (double x) override {
//        return -2 / (x*x*(x + 1));
//    }
//};
//
//class RealFunc: public Function1 {
//    double operator() (double x) {
//        return 1. / x + 1;
//    }
//};

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
    PFunc p;
    QFunc q;
    BVProblemL2R3 task(&f, left, right, lz0, re);
    
    ShootingMethodL2R3 s(task);
    s.Calc(h, eps);
    std::cout << '\n' << s.GetY() << '\n';

    FDMethod fd(&p, &q, task);
    fd.Calc(h);
    std::cout << '\n' << fd.GetY() << '\n';

    RealFunc rf;
    vdouble x = GetValInRange(left, right, h);
    vdouble ry = GetFValInRange(x, &rf);
    std::cout << ry;


    return 0;
}
