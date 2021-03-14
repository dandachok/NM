#include <iostream>
#include <vector>

#include "lib/NMLib.hpp"

double f(double x) {
    return x / std::pow(3*x + 4, 3);
}

int main() {
    double a = -1;
    double b = 1;
    double h1 = 0.5;
    double h2 = 0.25;

    Integrate intg(&f);

    double rech1 = intg.RectangleMethod(a, b, h1);
    double traph1 = intg.Trapezium(a, b, h1);
    double simph1 = intg.Simpson(a, b, h1);
    double rech2 = intg.RectangleMethod(a, b, h2);
    double traph2 = intg.Trapezium(a, b, h2);
    double simph2 = intg.Simpson(a, b, h2);
    double recrr = intg.RungeRombergOpt(rech1, rech2, h1, h2);
    double traprr = intg.RungeRombergOpt(traph1, traph2, h1, h2);
    double simprr = intg.RungeRombergOpt(simph1, simph2, h1, h2);

    std::cout << "\t\t\tH1\t\tH2\tRunge-Romberg\n";
    std::cout << "Rectangle method: " << rech1 << '\t' << rech2 << '\t' << recrr << '\n';
    std::cout << "Trapezium method: " << traph1 << '\t' << traph2 << '\t' << traprr << '\n';
    std::cout << "Simpson method: " << simph1 << '\t' << simph2 << '\t' << simprr << '\n';

    return 0;
    
}