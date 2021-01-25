#include "lib/NMLib.hpp"

int main () {
    int n;
    std::cin >> n;
    vvdouble mat(n, vdouble(n - 1));

    std::cin >> mat;

    vdouble ans = SweepMethod(mat); 
    std::cout << ans;

    return 0;
}