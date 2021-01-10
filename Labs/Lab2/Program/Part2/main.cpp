#include <iostream>
#include <vector>

using point = std::pair<double, double>;

double f1 (double x, double y) {
    return 0.1*x*x + x + 0.2*y*y - 0.3;
}

double f2 (double x, double y) {
    return 0.2*x*x + y - 0.1*x*y - 0.7;
}

double f1dx1 (double x, double y) {
    return 0.2*x + 1;
}

double f2dx1 (double x, double y) {
    return 0.4*x - 0.1*y;
}

double f1dx2 (double x, double y) {
    return 0.4*y;
}

double f2dx2 (double x, double y) {
    return 1 - 0.1*x;
}

double det (double a, double b, double c, double d) {
    return a * d - b * c;
}

double Norma (double nx, double ny, double x, double y) {
    return std::max(std::abs(nx - x), std::abs(ny - y));
}

double phi1 (double x, double y) {
    return 0.3 - 0.1*x*x - 0.2*y*y;
}

double phi2 (double x, double y) {
    return 0.7 - 0.2*x*x + 0.1*x*y;
}

point NewtonMethod (double x, double y, double eps) {
    double next_x = x;
    double next_y = y;
    do {
        x = next_x;
        y = next_y;
        double detJ = det(f1dx1(x,y), f1dx2(x,y), f2dx1(x,y), f2dx2(x,y));
        double detA_x = det(f1(x,y), f1dx2(x,y), f2(x,y), f2dx2(x,y));
        next_x = x - detA_x / detJ;

        double detA_y = det(f1dx1(x,y), f1(x,y), f1dx2(x,y), f2(x,y));
        next_y = y - detA_y / detJ;
    } while(Norma(next_x, next_y, x, y) > eps);

    return point(next_x, next_y);
}

point SimpleIter (double x, double y, double eps, double q) {
    double next_x = x;
    double next_y = y;
    do {
        x = next_x;
        y = next_y;
        next_x = phi1(x, y);
        next_y = phi2(x, y);
    } while (q / (1 - q) * Norma(next_x, next_y, x, y) > eps);
    return point(next_x, next_y);
}

int main() {
    double eps = 0.0001;
    double x = 0.25;
    double y = 0.75;
    double q = 0.5;
    //std::cout << (eps > 0);
    point ans = SimpleIter(x, y, eps, q);
    std::cout << ans.first << '\n';
    std::cout << ans.second << '\n';

    return 0;
}
