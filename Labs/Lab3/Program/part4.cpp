#include <iostream>
#include <vector>

#include "lib/NMLib.hpp"

int main() {
    
    vdouble x = {1.0, 1.5, 2.0, 2.5, 3.0};
    vdouble y = {0.0, 0.40547, 0.69315, 0.91629, 1.0986};
    double x_ = 2.0;

    Diff df(x, y);

    double d1 = df.dif(x_);
    double d2 = df.dif2(x_);

    std::cout << "First dif: " << d1 << '\n';
    std::cout << "Second dif " << d2 << '\n';

    return 0;
    
}