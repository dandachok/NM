#include <cmath>
#include <iostream>
#include <vector>

const double e = 2.718281828;
const double q = 0.64;

double f (double x) {
    return std::sqrt(1 - x*x) - std::pow(e, x) + 0.1;
}

double difF (double x) {
    return -x / std::sqrt(1 - x*x) - std::pow(e, x);
}

double NewtonMethod (double x0, double eps, double (*f)(double), double(*difF)(double)) {
    double x = x0;
    double next_x = x0;
    do {
        x = next_x;
        next_x = x - f(x)/difF(x);
    } while (std::abs(next_x - x) > eps);

    return next_x;
}

double phi (double x) {
    return std::log(std::sqrt(1 - x*x) + 0.1);
}

double SimpleIterMethod (double x0, double eps, double (*f)(double)) {
    double x = x0;
    double next_x = x0;
    do {
        x = next_x;
        next_x = f(x);
    } while (q / (1 - q) * std::abs(x - next_x) > eps);
    return next_x;
}

int main() {

    double Newton_x0 = 0.6;
    double SI_x0 = 0.475;
    double eps = 0.001;
    std::cout << "Newton method: " << NewtonMethod(Newton_x0, eps, f, difF) << '\n';
    std::cout << "Simple iter method: " << SimpleIterMethod(SI_x0, eps, phi) << '\n';

    return 0;
}