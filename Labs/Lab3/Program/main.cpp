#include <iostream>
#include <vector>

#include "TPolynomial.hpp"

double fLagrang (double x) {
    return log(x);
}

double fNewton (double x) {
    return sin(pi / 6 * x);
}

double fint (double x) {
    return x / pow(3*x + 4, 3);
}

int main() {
    int n;
    std::cin >> n;
    vdouble x(n);
    std::cin >> x;
    vdouble y(n);
    std::cin >> y;
    
    //double l = PolynomLagrange(x, f, x0);
    //double y = f(x0);
    //std::cout << "Lagrang polynomial res: " << l << '\n';
    //std::cout << "f(x) res: " << y << '\n';
    //std::cout << "diff: " << std::abs(l - y) << '\n';

    Integrate p(fint);
    double l = -1;
    double r = 1;
    double h1 = 0.25;
    double h2 = 0.5;
    double r1 = p.RectangleMethod(l, r, h1);
    double t1 = p.Trapezium(l, r, h1);
    double s1 = p.Simpson(l, r, h1);
    double r2 = p.RectangleMethod(l, r, h2);
    double t2 = p.Trapezium(l, r, h2);
    double s2 = p.Simpson(l, r, h2);
    std::cout << "H1\t\tH2\t\tRRR\n";
    std::cout << r1 << "\t" << r2 << '\t' << p.RungeRombergOpt(r1, r2, h1, h2) << '\n';
    std::cout << t1 << "\t" << t2 << '\t' << p.RungeRombergOpt(t1, t2, h1, h2) << '\n';
    std::cout << s1 << "\t" << s2 << '\t' << p.RungeRombergOpt(s1, s2, h1, h2) << '\n';

    //p.AddPoint(x0, fNewton(x0));
    //for (int i = 0; i < n; ++i) {
    //    std::cout << p(x[i]) << '\t' << y[i] << '\n';
    //}
    //std::cout << "SSE: " << p.SSE() << '\n';
    return 0;
}