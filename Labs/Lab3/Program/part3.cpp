#include <iostream>
#include <vector>

#include "lib/NMLib.hpp"

int main() {
    
    vdouble x = {-0.9, 0.0, 0.9, 1.8, 2.7, 3.6};
    vdouble y = {-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.313813};
    MNK f1(x, y, 1);
    MNK f2(x, y, 2);
    vdouble y1 = f1(x);
    vdouble y2 = f2(x);

    std::cout << "SSE first pow:" << f1.SSE() << '\n';
    std::cout << "SSE second pow: " << f2.SSE() << '\n';
    std::cout << "First pow result: " << y1 << '\n';
    std::cout << "Second pow result: " << y2 << '\n';
}