#include "lib/NMLib.hpp"

class MyFunc: public Function3 {
    double operator() (double x, double y, double z) override {
        return 2*y / (x*x*(x + 1));
    }
};

class PFunc: public Function1 {
    double operator() (double x) override {
        return 0;
    }
};

class QFunc: public Function1 {
    double operator() (double x) override {
        return -2 / (x*x*(x + 1));
    }
};

class RealFunc: public Function1 {
    double operator() (double x) {
        return 1. / x + 1;
    }
};

int main() {
    double left = 1.;
    double right = 2.;
    double z0 = -1;

    LinearDiffEqual de(2, -4, 4);
    MyFunc* f = new MyFunc();
    BVProblemL2R3 task(f, left, right, z0, de);
    
    ShootingMethodL2R3 sm(task);

    PFunc* p = new PFunc();
    QFunc* q = new QFunc();

    FDMethod fdm(p, q, task);

    double h = 0.01;
    double e = 0.001;

    auto x = GetValInRange(left, right, h);
    auto y = GetFValInRange(x, new RealFunc());

    std::cout << "Real function:\n";
    std::cout << y << '\n';

    sm.Calc(h, e);
    fdm.Calc(h);

    std::cout << "Shooting method:\n";
    std::cout << sm.GetY() << '\n';
    std::cout << "Finite-difference method:\n";
    std::cout << fdm.GetY() << '\n';

    return 0;
}