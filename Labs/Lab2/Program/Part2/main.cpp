#include <iostream>
#include <vector>
#include <cmath>

using point = std::pair<double, double>;

double f1 (double x, double y) {
    return (x*x + 16)*y - 64;
}

double f2 (double x, double y) {
    return (x - 2)*(x - 2) + (y - 2) * (y - 2) - 16;
}

double f1dx1 (double x, double y) {
    return 2*x*y;
}

double f2dx1 (double x, double y) {
    return 2*x - 4;
}

double f1dx2 (double x, double y) {
    return x*x + 16;
}

double f2dx2 (double x, double y) {
    return 2*y - 4;
}

double det (double a, double b, double c, double d) {
    return a * d - b * c;
}

double Norma (double nx, double ny, double x, double y) {
    return std::max(std::abs(nx - x), std::abs(ny - y));
}

double phi1 (double x, double y) {
    return std::sqrt(64/(x*x + 16));
}

double phi2 (double x, double y) {
    return std::sqrt(16 - (y - 2)*(y - 2)) + 2;
}

point NewtonMethod (double x, double y, double eps) {
    double next_x = x;
    double next_y = y;
    do {
        x = next_x;
        y = next_y;
        //std::cout << x << ' ' << y << '\n';
        double detJ = det(f1dx1(x,y), f1dx2(x,y), f1dx2(x,y), f2dx2(x,y));
        double detA_x = det(f1(x,y), f1dx2(x,y), f2(x,y), f2dx2(x,y));
        next_x = x - detA_x / detJ;

        double detA_y = det(f1dx1(x,y), f1(x,y), f2dx1(x,y), f2(x,y));
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
        next_x = phi2(x, y);
        next_y = phi1(x, y);
    } while (q / (1 - q) * Norma(next_x, next_y, x, y) > eps);
    return point(next_x, next_y);
}

int main() {
    double eps = 0.0001;
    double x = -1.7;
    double y = 3.3;
    double q = 0.6;
    //std::cout << (eps > 0);
    point ans = NewtonMethod(x, y, eps);
    std::cout << "Newton method: (" << ans.first << ", ";
    std::cout << ans.second << ")\n";
    ans = SimpleIter(x, y, eps, q);
    std::cout << "Simple iter method: (" << ans.first << ", ";
    std::cout << ans.second << ")\n";

    return 0;
}
