#include <iostream>
#include <vector>

#include "lib/NMLib.hpp"


int main() {
    
    vdouble x = {0.0, 0.9, 1.8, 2.7, 3.6};
    vdouble y = {0.0, 0.36892, 0.85408, 1.7856, 6.3138};
    double x0 = 1.5;

    CubicSplines cs(x, y);
    
    double yres = cs(1.5);

    std::cout << "Cubic spline in X* = " << x0 << ": " << yres << "\n";
    
    return 0;
}