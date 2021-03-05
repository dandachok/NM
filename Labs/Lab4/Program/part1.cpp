#include "lib/NMLib.hpp"

const double e = 2.7;

class MyFunc final: public Function3 {
    double operator() (double x, double y, double z) override {
        return 2*y + 4*x*x*std::pow(e, x*x);
    }
};

class RealFunc final: public Function1 {
    double operator() (double x) {
        return std::pow(e, x*x) + std::pow(e, x * std::sqrt(2)) + std::pow(e, -x * std::sqrt(2));
    }
};

int main() {
    double left = 0;
    double right = 1;
    double y0 = 3;
    double z0 = 0;
    MyFunc* f = new MyFunc();
    CauchyProblem cp(f, left, right, y0, z0);

    EilerMethod2 em(f, left, right, y0, z0);
    RungeKutta2 rg(cp);
    AdamsMethod am(cp);

    double h = 0.1;

    em(h);
    rg(h);
    am(h);

    std::cout << "Eiler method:\n";
    std::cout << em.GetY() << '\n';
    std::cout << "Runge-Kutta method:\n";
    std::cout << rg.GetY() << '\n';
    std::cout << "Adams method:\n";
    std::cout << am.GetY() << '\n';
    auto x = GetValInRange(left, right, h);
    auto y = GetFValInRange(x, new RealFunc());

    std::cout << "Real function:\n";
    std::cout << y << '\n';
    return 0;
}