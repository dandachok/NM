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
    am(h);
    auto em_y = em.GetY();
    auto rg_y1 = rg(h);
    auto rg_y2 = rg(h * 2);
    auto rr_y = RungeRomberg(rg_y1, rg_y2);
    auto am_y = am.GetY();

    auto x = GetValInRange(left, right, h);
    auto y = GetFValInRange(x, new RealFunc());

    std::cout << "Real\tEiler\tRK\tAdams\tRK(RR)\tDelta\n";
    for (int i = 0; i < em_y.size(); ++i) {
        std::cout << y[i] << '\t' << em_y[i] << '\t' << rg_y1[i] << '\t';
        std::cout << am_y[i];
        if (i % 2 == 0) std::cout << '\t' << rr_y[i/2] << '\t' << std::abs(y[i] - rr_y[i/2]);
        std::cout << '\n';
    }

    //std::cout << "Eiler method:\n";
    //std::cout << em.GetY() << '\n';
    //std::cout << "Runge-Kutta method:\n";
    //std::cout << rg_y1 << '\n' << "RR optimization:";
    //std::cout << RungeRomberg(rg_y1, rg_y2) << '\n';
    //rg(h*2);
    //std::cout << "Adams method:\n";
    //std::cout << am.GetY() << '\n';

    return 0;
}