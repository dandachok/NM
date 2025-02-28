#include "lib/NMLib.hpp"

int main () {
    int n;
    std::cin >> n;
    vvdouble mat(n, vdouble(n));

    std::cin >> mat;
    vdouble b(n);
    std::cin >> b;
    vvdouble ludec = LUDecomposition(mat);
    vvdouble l = GetL(ludec);
    vvdouble u = GetU(ludec);

    double det = Det(ludec);

    vvdouble inv = InversMatrix(ludec);

    std::cout << "L: " << l << "U:" << u;
    vdouble ans = LUSolve(ludec, b);
    std::cout << "LU solve:" << ans << '\n';
    ans = GausseMethod(mat, b);
    std::cout << "Gausse solve:" << ans << '\n';
    std::cout << "Determinant: " << det << '\n';
    std::cout << "Invert matrix:\n" << inv << '\n';
    std::cout << mat * inv << '\n';
    return 0;
}
