#include "lib/NMLib.hpp"

int main () {
    int n;
    std::cin >> n;

    vvdouble mat(n, vdouble(n));
    std::cin >> mat;
    double eps = 0.01;
    vvdouble ans = QRMethod(mat, eps);
    std::cout << ans;
    return 0;
}