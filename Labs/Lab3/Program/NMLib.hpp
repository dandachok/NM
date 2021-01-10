#pragma once

#include <cmath>
#include <iostream>
#include <vector>

using complex = std::pair<double, double>;
using vcomplex = std::vector<complex>;
using vint = std::vector<int>;
using vvint = std::vector<std::vector<int>>;
using vdouble = std::vector<double>;
using vvdouble = std::vector<std::vector<double>>;

const double pi = 3.1415926535;

double sign (double v) {
    if (v == 0) {
        return 0;
    }
    return v > 0? 1: -1;
}

// Matrix operations

std::istream& operator>> (std::istream& in, vvdouble& mat) {
    for (auto& i: mat) {
        for (auto& j: i) {
            in >> j;
        }
    }
    return in;
}

std::istream& operator>> (std::istream& in, vdouble& v) {
    for (auto& i: v) {
        in >> i;
    }

    return in;
}

std::ostream& operator<< (std::ostream& out, const vvdouble& mat) {
    out.precision(3);
    for (const auto& i: mat) {
        for (const auto& j: i) {
            out <<  j << '\t';
        }

        out << '\n';
    }

    return out;
}

std::ostream& operator<< (std::ostream& out, const vdouble& mat) {
    for (const auto& i: mat) {
            out << i << ' ';
    }
    out << '\n';

    return out;
}

complex operator+ (complex a, const complex& c) {
    a.first += c.first;
    a.second += c.second;
    return a;
}

vdouble operator+ (vdouble a, const vdouble& b) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        a[i] += b[i];
    }
    return a;
}

vvdouble operator+ (vvdouble a, const vvdouble& b) {
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[i][j] += b[i][j];
        }
    }
    return a;
}

complex operator* (double a, complex b) {
    b.first *= a;
    b.second *= a;
    return b;
}

vvdouble operator* (double a, vvdouble b) {
    for (auto& i: b) {
        for (auto& j: i) {
            j *= a;
        }
    }
    return b;
}

vdouble operator* (vdouble a, const vdouble& b) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        a[i] *= b[i];
    }
    return a;
}

vdouble operator* (const vvdouble& a, const vdouble& b) {
    int n = a.size();
    int m = b.size();
    vdouble res(a.size(), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i] += a[i][j] * b[j];
        }
    }
    return res;
}

vvdouble operator* (const vvdouble& a, const vvdouble& b) {
    int n = a.size();
    int m = b[0].size();
    int l = a[0].size();

    vvdouble res(n, vdouble(m, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < l; ++k) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return res;
}

complex operator- (complex a, const complex& b) {
    return a + (-1) * b;
}

vdouble operator- (vdouble a, const vdouble& b) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        a[i] -= b[i];
    }
    return a;
}

vvdouble operator- (vvdouble a, const vvdouble& b) {
    return a + (-1)*b;
}

vdouble Pow (vdouble a, int k) {
    for (int i = 0; i < a.size(); ++i) {
        a[i] = pow(a[i], k);
    }
    return a;
}

double Sum (const vdouble& v) {
    double sum = 0;
    for (auto i: v) {
        sum += i;
    }
    return sum;
}

double Norma (const complex& c) {
    return sqrt(c.first * c.first + c.second * c.second);
}

double Norma (const vdouble& v) {
    int n = v.size();
    double norma = 0;
    for (int i = 0; i < n; ++i) {
        norma += std::abs(v[i]);
    }

    return norma;
}

double Norma (const vvdouble& mat) {
    int n = mat.size();
    int m = mat[0].size();
    double max = 0;
    for (int j = 0; j < m; ++j) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += std::abs(mat[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }

    return max;
}

double Norma2 (const vvdouble& mat) {
    double sum = 0;
    for (const auto& i: mat) {
        for (const auto& j: i) {
            sum += j * j;
        }
    }

    return sqrt(sum);
}

vvdouble Trans (const vvdouble& mat) {
    int m = mat.size();
    int n = mat[0].size();

    vvdouble res(n, vdouble(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i][j] = mat[j][i];
        }
    }
    return res;
}

vvdouble CreateIdentity (int n) {
    vvdouble im(n, vdouble(n, 0));
    for (int i = 0; i < n; ++i) {
        im[i][i] = 1;
    }

    return im;
}

// Algorithms

vdouble GausseMethod (vvdouble mat) {
    int n = mat.size();
    int m = mat[0].size() - 1;

    vdouble ans(n, 0);
    //forward
    for (int j = 0; j < std::min(n, m); ++j) {
        for (int i = j + 1; i < n; ++i) {
            double r = -mat[i][j] / mat[j][j];
            for (int k = j; k < m + 1; ++k) {
                mat[i][k] += mat[j][k] * r;
            }
        }
    }
    //back
    for (int i = n - 1; i >= 0; --i) {
        ans[i] = mat[i][m];
        for (int j = i + 1; j < n; ++j) {
            ans[i] -= ans[j] * mat[i][j];
        }
        ans[i] /= mat[i][i];
    }
    return ans;
}


vdouble SweepMethod (vvdouble mat) {
    int n = mat.size();
    vdouble p(n, 0);
    vdouble q(n, 0);
    p[0] = -mat[0][2] / mat[0][1];
    q[0] = mat[0][3] / mat[0][1];
    for (int i = 1; i < n; ++i) {
        double a = mat[i][0];
        double b = mat[i][1];
        double c = mat[i][2];
        double d = mat[i][3];
        p[i] = - c / (b + a*p[i - 1]);
        q[i] = (d - a*q[i-1]) / (b + a*p[i-1]);
    }
    p[n - 1] = 0;
    vdouble x(n);
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i]*x[i + 1] + q[i];
    }
    //std::cout <<"P: \n" << p << "Q\n" << q;
    return x;
}


vvdouble LUDecomposition (vvdouble mat) {
    int n = mat.size();
    int m = mat[0].size() - 1;

    vdouble ans(n, 0);
    vint swapMatrix(n);
    for (int j = 0; j < std::min(n, m); ++j) {
        
        double max = 0;
        int str;
        for (int i = j; i < n; ++i) {
            if (mat[i][j] > max) {
                max = mat[i][j];
                str = i;
            }
        }
        swapMatrix[j] = str;
        std::swap(mat[j], mat[str]);
        for (int i = j + 1; i < n; ++i) {
            double r = -mat[i][j] / mat[j][j];
            for (int k = j; k < m; ++k) {
                mat[i][k] += mat[j][k] * r;
            }
            mat[i][j] = -r;
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        std::swap(mat[i], mat[swapMatrix[i]]);
    }
    return mat;
}


vdouble SimpleIter (const vvdouble& m, vdouble b, double eps) {
    int n = m.size();
    vvdouble a = m;
    for (int i = 0; i < n; ++i) {
        b[i] /= m[i][i];
        for (int j = 0; j < n; ++j) {
            a[i][j] = i == j ? 0: -m[i][j] / m[i][i]; 
        }
    }

    vdouble x = b;
    vdouble prev_x;
    int count_iter = 0;
    do {
        count_iter++;
        prev_x = x;
        x = b + a * prev_x;
    } while (Norma(x - prev_x) > eps);

    //std::cout << count_iter << '\n';
    return x;
}


vdouble ZeidelMethod (const vvdouble& m, vdouble b, double eps) {
    int n = m.size();
    vvdouble a = m;
    for (int i = 0; i < n; ++i) {
        b[i] /= m[i][i];
        for (int j = 0; j < n; ++j) {
            a[i][j] = i == j ? 0: -m[i][j] / m[i][i]; 
        }
    }

    vdouble x = b;
    vdouble prev_x;
    double count_iter = 0;
    do {
        count_iter++;
        prev_x = x;
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j) {
                x[i] += x[j] * a[i][j];
            }
        }
    } while (Norma(x - prev_x) > eps);
    
    std::cout << count_iter << '\n';

    return x;
}


void FindMaxNDElem (const vvdouble& a, int& res_i, int& res_j) {
    int n = a.size();
    double max = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (a[i][j] > max && i != j) {
                max = a[i][j];
                res_i = i;
                res_j = j;
            }
        }
    }
}

double SqrtSumNDElem (const vvdouble& a) {
    int n = a.size();
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i < j) {
                sum += a[i][j] * a[i][j];
            }
        }
    }
    return sqrt(sum);
}

vvdouble JakobiMethod (const vvdouble& mat, double eps) {
    int n = mat.size();
    int im;
    int jm;
    vvdouble a = mat;
    for (int i = 0; SqrtSumNDElem(a) > eps; ++i) {

        FindMaxNDElem(a, im, jm); //find max not diagonal elem

        double q;
        if (a[im][im] == a[jm][jm]) {
            q = pi / 4; 
        } else {
            q = 1./2 * atan (2*a[im][jm] / (a[im][im] - a[jm][jm]));
        }

        //if (i < 10) {
        //    std::cout << im << ' ' << jm << ' ' << q << '\n';
//
        //}

        vvdouble u = CreateIdentity(n);
        u[im][im] = cos(q);
        u[im][jm] = -sin(q);
        u[jm][im] = sin(q);
        u[jm][jm] = cos(q);
        vvdouble ut = Trans(u);

        //std::cout << "U:\n" << u;
        
        a = ut * a * u;
    }

    return a;
}


vvdouble HausseholderMatrix (const vvdouble& v) {
    return CreateIdentity(v.size()) - 2 / (Trans(v) * v)[0][0] * (v * Trans(v));
}

vvdouble GetColumn (const vvdouble& mat, int j) {
    int n = mat.size();
    double norma = 0;
    for (int i = j; i < n; ++i) {
        norma += mat[i][j] * mat[i][j];
    }
    norma = sqrt(norma);

    vvdouble clm(n, vdouble(1, 0));
    
    clm[j][0] = mat[j][j] + sign(mat[j][j]) * norma;
    for (int i = j + 1; i < n; ++i) {
        clm[i][0] = mat[i][j];
    }
    return clm;
}

void QRDecomposition (vvdouble a, vvdouble& q, vvdouble& r) {
    int n = a.size();
    q = CreateIdentity(n);
    for (int i = 0; i < n - 1; ++i) {
        vvdouble b = GetColumn(a, i);
        vvdouble h = HausseholderMatrix(b);
        a = h * a;
        q = q * h;
    }
    r = a;
}

double SqrtSumSubStr (const vvdouble& a, int i, int j) {
    int n = a.size();
    double sum = 0;
    for (int k = i + 1; k < n; ++k) {
        sum += a[k][j] * a[k][j];
    }
    return sqrt(sum);
}

complex Solve (double a, double b, double c) {
    double d = b*b - 4*a*c;
    complex ans;
    if (d < 0) {
        ans.first = -b / 2 / a;
        ans.second = sqrt(abs(d)) / 2 / a;
    } else {
        ans.first = (-b + sqrt(d)) / 2 / a;
        ans.second = 0;
    }
    return ans;
}

bool FinishIterProc(const vvdouble& a, const vvdouble& ap, double e) {
    int n = a.size();
    for (int i = 0; i < n - 1; ++i) {
        if (SqrtSumSubStr(a, i, i) > e) {
            if (SqrtSumSubStr(a, i + 1, i) > e) {
                return false;
            } else {
                double b = -a[i][i]-a[i+1][i+1];
                double c = -a[i+1][i] * a[i][i+1];
                complex l = Solve (1, b, c);
                b = -ap[i][i]-ap[i+1][i+1];
                c = -ap[i+1][i] * ap[i][i+1];
                complex lp = Solve(1, b, c);
                if (Norma(l - lp) > e) {
                    return false;
                }
            }
        }
    }
    return true;
}

vvdouble QRMethod (vvdouble a, double eps) {
    int n = a.size();
    vvdouble q, r, an = a;
    vcomplex l(n);
    int i = 0;
    do {
        std::swap(a, an);
        QRDecomposition(a, q, r);
        an = r * q;
    } while (!FinishIterProc(an, a, eps));
    return an;
}
