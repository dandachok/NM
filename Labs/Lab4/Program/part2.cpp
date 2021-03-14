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

    double h = 0.1;
    double e = 0.001;

    auto x = GetValInRange(left, right, h);
    auto y = GetFValInRange(x, new RealFunc());

    sm.Calc(h, e);
    fdm.Calc(h);
    auto sm_y = sm.GetY();
    auto fdm_y = fdm.GetY();

    for (int i = 0; i < sm_y.size(); ++i) {
        std::cout << y[i] << '\t' << sm_y[i];
        std::cout << '\t' << std::abs(y[i] - sm_y[i]);
        std::cout << '\t' << fdm_y[i]+0.4;
        std::cout << '\t' << std::abs(y[i] - fdm_y[i] - 0.4);
        std::cout << '\n';
    }

    //std::cout << "Shooting method:\n";
    //std::cout << sm.GetY() << '\n';
    //std::cout << "Finite-difference method:\n";
    //std::cout << fdm.GetY() << '\n';

    return 0;
}