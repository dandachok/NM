#include <iostream>
#include <vector>

#include "lib/NMLib.hpp"

double f (double x) {
    return std::tan(x);
}

double fNewton (double x) {
    return sin(pi / 6 * x);
}

double fint (double x) {
    return x / pow(3*x + 4, 3);
}

int main() {
    
    vdouble x = {0., pi/8, 2*pi/8, 3*pi/8};
    double x0 = 3*pi/16;
    double real = f(x0);
    std::cout << "X*:\t" << x0 << '\n';
    std::cout << "Real:\t" << real << '\n';
    std::cout << '\n';

    std::cout << "Test case 1:\t" << x << '\n';
    double plr1 = PolynomLagrange(x, f, x0);
    std::cout << "Lagrange method: " << plr1 << ", Error: ";
    std::cout << std::abs(plr1 - real) << '\n';

    PolynomNewton pn(x, GetY(x, f));
    double pnr1 = pn(x0);
    std::cout << "Newton method: " << pnr1 << ", Error: ";
    std::cout << std::abs(pnr1 - real) << '\n';

    std::cout << '\n';

    x = {0., pi/8, pi/3, 3*pi/8};
    x0 = 3*pi/16;
    std::cout << "Test case 2:\t" << x << '\n';
    double plr2 = PolynomLagrange(x, f, x0);
    std::cout << "Lagrange method: " << plr2 << ", Error: ";
    std::cout << std::abs(plr2 - real) << '\n';
    PolynomNewton pn2(x, GetY(x, f));
    double pnr2 = pn2(x0);
    std::cout << "Newton method: " << pnr2 << ", Error: ";
    std::cout << std::abs(pnr2 - real) << '\n';
    
    
    return 0;
}