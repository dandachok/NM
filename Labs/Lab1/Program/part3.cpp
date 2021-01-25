#include "lib/NMLib.hpp"

int main () {
    int n, m;

    std::cin >> n >> m;

    vvdouble A(n, vdouble(m));
    vdouble b(n);

    std::cin >> A >> b;

    double eps;
    std::cin >> eps;

    vdouble simp_ans = SimpleIter(A, b, eps);
    vdouble zeidel_ans = ZeidelMethod(A, b, eps);
    std::cout << "Simple ans:\n" << simp_ans;
    std::cout << "Zeidel ans:\n" << zeidel_ans;
    return 0;
}