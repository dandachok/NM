#include "lib/NMLib.hpp"

int main () {
    int n;
    std::cin >> n;

    vvdouble mat(n, vdouble(n));
    std::cin >> mat;
    double eps = 0.01;
    vdouble ans = JakobiMethod(mat, eps);
    std::cout << ans;
    return 0;
}