#pragma once

#include "FunAlgo.hpp"

using dbl = double;

class Function1 {
    public:
    virtual double operator() (double) = 0;
};

class Function3 {
    public:
    virtual double operator() (double, double, double) = 0;

    virtual double c (double a, double b, double c) {
        return operator() (a, b, c);
    }
};

vdouble RungeRomberg(vdouble y1, vdouble y2) {
    vdouble ans;
    int n = y1.size();
    for (int i = 0; i < n; i += 2) {
        double res = y1[i] + (y1[i] - y2[i/2])/15;
        ans.push_back(res);
    }
    return ans;
}

// LDE has the following pattern: a*f'(x) + b*f(x) = c

class LinearDiffEqual {
    public:
        LinearDiffEqual() = default;
        LinearDiffEqual(const LinearDiffEqual&) = default;
        LinearDiffEqual(double a, double b, double c) :
            a(a), b(b), c(c) {}
        
        LinearDiffEqual& operator= (const LinearDiffEqual&) = default;

        double a;
        double b;
        double c;
};

class CauchyProblem {
    public:
        CauchyProblem() = default;
        CauchyProblem(const CauchyProblem& de) = default;
        CauchyProblem(Function3* f, double left, double right, double y0, double y1) :
            f(f), y0(y0), z0(y1), left(left), right(right) {}

        CauchyProblem& operator= (const CauchyProblem&) = default;

        Function3* f;
        double y0;
        double z0;
        double left;
        double right;
};

using LDE = LinearDiffEqual;
using CP = CauchyProblem;

class BVProblemL2R3: public CauchyProblem {
    public:
        BVProblemL2R3() = default;
        BVProblemL2R3(const BVProblemL2R3&) = default;
        BVProblemL2R3(Function3* f, double left, double right,
                    double lz0, LDE e2) :
            CauchyProblem(f, left, right, 0, lz0), re(e2) {}
        
        BVProblemL2R3& operator= (const BVProblemL2R3&) = default;

        LDE re;
};

vdouble GetValInRange(double left, double right, double h) {
    vdouble res;
    double eps = 0.0000001;
    while (left < right || (std::abs(right - left) < eps)) {
        res.push_back(left);
        left += h;
    }
    return res;
}

vdouble GetFValInRange(const vdouble& x, Function1* f) {
    vdouble y;
    for (int i = 0; i < x.size(); ++i) {
        y.push_back((*f)(x[i]));
    }
    return y;
}

class EilerMethod2 {
    public:
        EilerMethod2(Function3* f, double left, double right, double y0, double z0) :
            left(left),
            right(right),
            f(f),
            y0(y0),
            z0(z0) {
        }
    
        void operator() (double h) {
            x = GetValInRange(left, right, h);
            y.clear();
            z.clear();
            y.push_back(y0);
            z.push_back(z0);

            int n = x.size() - 1;
            for (int i = 0; i < n; ++i) {
                double z_next = z[i] + h * (*f)(x[i], y[i], z[i]);
                z.push_back(z_next);
                double y_next = y[i] + h * z[i];
                y.push_back(y_next);
            }
        }

        vdouble GetY() {
            return y;
        }

        vdouble GetZ() {
            return z;
        }
    

    private:
        double left;
        double right;
        Function3* f;
        double h;
        double y0;
        double z0;

        vdouble x;
        vdouble y;
        vdouble z;
};

class RungeKutta2 {
    public:
        RungeKutta2(CauchyProblem t) :
            left(t.left),
            right(t.right),
            f(t.f),
            y0(t.y0),
            z0(t.z0) {
        }
        
        vdouble operator() (double h) {
            x = GetValInRange(left, right, h);
            int n = x.size();
            y.clear();
            z.clear();
            y.push_back(y0);
            z.push_back(z0);

            for (int i = 0; i < n - 1; ++i) {
                vdouble k(ORDER);
                vdouble l(ORDER);
                for (int j = 0; j < ORDER; ++j) {
                    if (j == 0) {
                        l[j] = h * (*f)(x[i], y[i], z[i]);
                        k[j] = h * z[i];
                    }
                    else if (j == 3) {
                        l[j] = h * (*f)(x[i] + h, y[i] + k[j - 1], z[i] + l[j - 1]);
                        k[j] = h * (z[i] + l[j - 1]);
                    } else {
                        l[j] = h * (*f)(x[i] + 0.5*h, y[i] + 0.5*k[j - 1], z[i] + 0.5*l[j - 1]);
                        k[j] = h * (z[i] + 0.5*l[j - 1]);
                    }
                }

                double y_next = y[i] + dif(k);
                double z_next = z[i] + dif(l);
                y.push_back(y_next);
                z.push_back(z_next);
            }

            return y;
        }

        vdouble GetY() {
            return y;
        }

        vdouble GetZ() {
            return z;
        }

        vdouble Get() {
            return y;
        }
        
        Function3* f;
        double left;
        double right;
        double y0;
        double z0;

    private:
        double dif (const vdouble& k) {
            return (k[0] + 2*(k[1] + k[2]) + k[3]) / 6;
        }

        double g (double z) {
            return z;
        }

        static const int ORDER = 4;
        
        vdouble x;
        vdouble y;
        vdouble z;
};

class AdamsMethod {
    public:
        AdamsMethod(CauchyProblem t) :
            f(t.f),
            left(t.left),
            right(t.right),
            y0(t.y0),
            z0(t.z0),
            task(t) {}
        
        void operator() (double h) {
            RungeKutta2 rk(task);
            rk(h);
            x = GetValInRange(left, right, h);
            y = rk.GetY();
            z = rk.GetZ();

            y.resize(4);
            z.resize(4);

            int n = x.size();
            for (int i = 3; i < n - 1; ++i) {
                double z_next = z[i] + h / 24 * df(i);
                z.push_back(z_next);
                double y_next = y[i] + h / 24 * dg(i);
                y.push_back(y_next);
            }
        }

        vdouble GetY() {
            return y;
        }

        vdouble GetZ() {
            return z;
        }
    private:

        double df(int i) {
            double f0 = (*f)(x[i], y[i], z[i]);
            double f1 = (*f)(x[i - 1], y[i - 1], z[i - 1]);
            double f2 = (*f)(x[i - 2], y[i - 2], z[i - 2]);
            double f3 = (*f)(x[i - 3], y[i - 3], z[i - 3]);

            return (55*f0 - 59*f1 + 37*f2 - 9*f3);
        }

        double dg(int i) {
            return (55*z[i] - 59*z[i - 1] + 37*z[i - 2] - 9*z[i - 3]);
        }

        CauchyProblem task;
        Function3* f;
        double left;
        double right;
        double y0;
        double z0;

        vdouble x;
        vdouble y;
        vdouble z;
};

class ShootingMethodL2R3 {
    public:
        ShootingMethodL2R3(BVProblemL2R3 task) : task(task), f(task.f) {}

        void Calc (double h, double eps) {
            x = GetValInRange(task.left, task.right, h);

            CauchyProblem ctask1 = task;
            CauchyProblem ctask2 = task;

            ctask1.y0 = 100;
            ctask2.y0 = 1000;

            RungeKutta2 rg1(ctask1);
            RungeKutta2 rg2(ctask2);

            double& n1 = rg1.y0;
            double& n2 = rg2.y0;

            vdouble y1 = rg1(h);
            vdouble y2 = rg2(h);

            if (isFinish(y1, eps)) {
                y = y1;
                return;
            } else if (isFinish(y2, eps)) {
                y = y2;
                return;
            }
            do {
                

                double next_n = GetNextN(n1, n2, y1, y2);
                n1 = n2;
                n2 = next_n;
                y1 = y2;
                y2 = rg2(h);
            } while (!isFinish(y2, eps));

            y = y2;
        }

        vdouble GetY() const {
            return y;
        }

    private:

        double GetYDer(const vdouble& y, double xp) {
            int i = 0;
            for (; i < x.size() - 2; ++i) {
                if (x[i] <= xp && xp <= x[i + 1]) {
                    break;
                }
            }
            return (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }

        double GetPhi (const vdouble& y) {
            double y_der = GetYDer(y, task.right);
            double phi = task.re.b * y_der + task.re.a * y[y.size() - 1] - task.re.c;
            return phi;
        }

        double GetNextN(double n1, double n2, const vdouble& y1, const vdouble& y2) {
            double phi1 = GetPhi(y1);
            double phi2 = GetPhi(y2);

            return n2 - (n2 - n1) / (phi2 - phi1) * phi2;
        }

        bool isFinish(const vdouble& y, double eps) {
            double phi = GetPhi(y);
            return std::abs(phi) < eps;
        }

        Function3* f;
        BVProblemL2R3 task;

        vdouble x;
        vdouble y;
};

class FDMethod {
    public:
    FDMethod(Function1* p, Function1* q, BVProblemL2R3 t) :
        p(p), q(q), task(t) {}

    void Calc(double h) {
        x = GetValInRange(task.left, task.right, h);
        int n = (task.right - task.left) / h;
        vvdouble mat(n + 1, vdouble(4));
        mat[0][0] = 0;
        mat[0][1] = -1;
        mat[0][2] = 1;
        mat[0][3] = task.z0 * h;


        for (int i = 1; i < n; ++i) {
            mat[i][0] = 1 - (*p)(x[i]) * h / 2;
            mat[i][1] = (*q)(x[i]) * h * h - 2;
            mat[i][2] = 1 + (*p)(x[i]) * h / 2;
            mat[i][3] = 0;
        }


        mat[n][0] = -task.re.b;
        mat[n][1] = h * task.re.a + task.re.b;
        mat[n][2] = 0;
        mat[n][3] = task.re.c * h;


        y = SweepMethod(mat);
    }

    vdouble GetY() {
        return y;
    }


    private:
        Function1* p;
        Function1* q;
        BVProblemL2R3 task;

        vdouble x;
        vdouble y;
};
