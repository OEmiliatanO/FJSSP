#include <algorithm>

class TFN_t
{
public:
    double l = 0, m = 0, r = 0;

    TFN_t() = default;
    TFN_t(double l_, double m_, double r_): l(l_), m(m_), r(r_) {};
    TFN_t(const TFN_t& other): l(other.l), m(other.m), r(other.r) {};

private:
};

TFN_t operator+(const TFN_t& A, const TFN_t& B)
{
    return TFN_t{A.l + B.l, A.m + B.m, A.r + B.r};
}

TFN_t max(const TFN_t& A, const TFN_t& B)
{
    return TFN_t{std::max(A.l, B.l), std::max(A.m, B.m), std::max(A.r, B.r)};
}

double PSD(const TFN_t& B, const TFN_t& A)
{
    if (B.r <= A.m) return 0;
    if (B.m < A.r and B.r > A.m)
    {
        double x0 = (1.0 + A.r / (A.m - A.r) + B.r / (B.m - B.r)) / (1.0 / (A.m - A.r) + 1.0 / (B.m - B.r));
        return (x0 - B.r) / (B.m - B.r);
    }
    return 1;
}

double ND(const TFN_t& B, const TFN_t& A)
{
    if (A.l >= B.m) return 0;
    if (A.l < B.m and A.m > B.l)
    {
        double x0 = (1.0 + B.l / (B.m - B.l) + A.l / (A.m - A.l)) / (1.0 / (A.m - A.l) + 1.0 / (B.m - B.l));
        return (x0 - A.l) / (A.m - A.l);
    }
    return 1;
}
