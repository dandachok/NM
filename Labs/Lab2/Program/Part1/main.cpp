#include <cmath>
#include <iostream>
#include <vector>

const double e = 2.7;
const double q = 0.64;

double f (double x) {
    return std::pow(e, 2*x) + 3*x - 4;
}

double difF (double x) {
    return 2 * std::pow(e, 2*x) + 3;
}

double NewtonMethod (double x0, double e, double (*f)(double), double(*difF)(double)) {
    double x = x0;
    double next_x = x0;
    do {
        x = next_x;
        next_x = x - f(x)/difF(x);
    } while (std::abs(next_x - x) > e);

    return next_x;
}

double phi (double x) {
    return std::log(4 - 3*x) / 2;
}

double SimpleIterMethod (double x0, double e, double (*f)(double)) {
    double x = x0;
    double next_x = x0;
    do {
        x = next_x;
        next_x = f(x);
    } while (q / (1 - q) * std::abs(x - next_x) > e);
    return next_x;
}

int main() {

    double Newton_x0 = 0.6;
    double SI_x0 = 0.475;
    double e = 0.001;
    std::cout << NewtonMethod(Newton_x0, e, f, difF) << '\n';
    std::cout << SimpleIterMethod(SI_x0, e, phi) << '\n';

    return 0;
}