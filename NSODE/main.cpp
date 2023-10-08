#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <random>
#include "fuzzy.h"
#define SEED 100

constexpr int MAXN = 100, MAXM = 100;
constexpr int MAX_POP = 1000;

TFN_t P[MAXN][MAXM][MAXM], D[MAXN], C[MAXN];
int n, m;

int seed = SEED;
int popn, itn;
double CR, MC, omega1 = 0.5, omega2 = 0.5;

std::default_random_engine gen(time(NULL));
std::uniform_real_distribution<double> rnd01(0, 1);

std::vector<double> x[MAX_POP], d[MAX_POP], v[MAX_POP], best_sol;
double best_fitness = std::numeric_limits<double>::min();
std::vector<int> op_on_machine[MAXN];

void init()
{
    best_sol.resize(n * m + 1);
    for (int i = 1; i <= popn; ++i)
    {
        x[i].resize(n * m + 1);
        d[i].resize(n * m + 1);
        v[i].resize(n * m + 1);
        for (int j = 1; j <= n * m; ++j)
            x[i][j] = rnd01(gen);
    }
}

void mutation()
{
    for (int i = 1; i <= popn; ++i)
    {
        //d[i] = x[m1] + C(x[m2] - x[m3]);
        int m1 = static_cast<int>(rnd01(gen) * popn) % popn + 1, m2 = static_cast<int>(rnd01(gen) * popn) % popn + 1, m3 = static_cast<int>(rnd01(gen) * popn) % popn + 1;
        for (int j = 1; j <= n * m; ++j)
        {
            d[i][j] = x[m1][j] + MC * (x[m2][j] - x[m3][j]);
        }
    }
}

void crossover()
{
    for (int i = 1; i <= popn; ++i)
    {
        int randj = static_cast<int>(rnd01(gen) * popn);
        for (int j = 1; j <= n * m; ++j)
            v[i][j] = (rnd01(gen) <= CR or j == randj ? d[i][j] : x[i][j]); 
    }
}

void print_sol(const std::vector<double>& x)
{
    std::vector<std::pair<double, int>> tmp;
    tmp.resize(n * m + 1);
    int op[MAXM];
    std::fill(op, op + n + 1, 1);
    for (int i = 1; i <= n * m; ++i)
        tmp[i] = {x[i], i};
    TFN_t machine_t[MAXM]{}, cp[MAXN]{};
    std::sort(tmp.begin() + 1, tmp.end(), std::greater<std::pair<double, int>>());
    for (int i = 1; i <= n * m; ++i)
    {
        int j = tmp[i].second % n + 1;
        std::cout << j << ' ';
        int machine = op_on_machine[j][op[j]];
        cp[j] = max(cp[j], machine_t[machine]) + P[j][op[j]][machine];
        machine_t[machine] = cp[j];
        op[j]++;
    }
    std::cout << '\n';
    for (int i = 1; i <= n; ++i)
    {
        std::cout << "job " << i << ": PSD = " << PSD(D[i], cp[i]) << " " << ND(D[i], cp[i]) << '\n';
        std::cout << "job " << i << " = (" << cp[i].l << ' ' << cp[i].m << ' ' << cp[i].r << ")" << '\n';
        std::cout << "job " << i << " = (" << D[i].l << ' ' << D[i].m << ' ' << D[i].r << ")" << '\n';
    }
    std::cout << '\n';
}

double f(const std::vector<double>& x)
{
    std::vector<std::pair<double, int>> tmp;
    int op[MAXM];
    std::fill(op, op + n + 1, 1);
    tmp.resize(n * m + 1);
    for (int i = 1; i <= n * m; ++i)
        tmp[i] = {x[i], i};
    std::sort(tmp.begin() + 1, tmp.end(), std::greater<std::pair<double, int>>());
    /*for (auto& p : tmp)
        std::cerr << '(' << p.first << ", " << p.second << ") ";
    std::cerr << '\n';*/
    TFN_t machine_t[MAXM]{}, cp[MAXN]{};
    for (int i = 1; i <= n * m; ++i)
    {
        int j = tmp[i].second % n + 1;
        int machine = op_on_machine[j][op[j]];
        cp[j] = max(cp[j], machine_t[machine]) + P[j][op[j]][machine];
        machine_t[machine] = cp[j];
        op[j]++;
    }

    double obj = 0;
    for (int i = 1; i <= n; ++i)
    {
        //std::cout << PSD(D[i], cp[i]) << " " << ND(D[i], cp[i]) << '\n';
        obj += omega1 * PSD(D[i], cp[i]) + omega2 * ND(D[i], cp[i]);
    }
    
    return obj;
}

void selection()
{
    for (int i = 1; i <= popn; ++i)
    {
        double fvi = f(v[i]), fxi = f(x[i]);
        //std::cerr << fvi << ' ' << fxi << '\n';
        if (fvi > fxi)
        {
            std::copy(v[i].begin(), v[i].end(), x[i].begin());
            if (best_fitness < fvi)
            {
                best_sol = v[i];
                best_fitness = fvi;
            }
        }
        else
        {
            if (best_fitness < fxi)
            {
                best_sol = x[i];
                best_fitness = fxi;
            }
        }
    }
}

void improved_selection()
{
    std::pair<double, int> tmp[MAXN * 2 + 1];
    std::vector<double> S[MAXN * 2 + 1];
    for (int i = 1; i <= popn; ++i)
    {
        S[i] = x[i];
        S[popn + i] = v[i];
        tmp[i] = {f(x[i]), i};
        tmp[popn + i] = std::make_pair<double, int>(f(v[i]), popn + i);
    }
    std::sort(tmp + 1, tmp + popn * 2 + 1);
    if (best_fitness < tmp[popn * 2].first)
    {
        best_sol = S[tmp[popn * 2].second];
        best_fitness = tmp[popn * 2].first;
    }
    for (int i = popn * 2; i >= popn + 1; --i)
        x[i] = S[tmp[i].second];
}

int main()
{
    std::cin >> popn >> itn >> CR >> MC;
    std::cout << popn << ' ' << itn << ' ' << CR << ' ' << MC << '\n';

    std::cin >> n >> m;

    for (int i = 1; i <= n; ++i)
    {
        op_on_machine[i].emplace_back(-1); // dummy
        for (int j = 1, machine; j <= m; ++j)
        {
            std::cin >> machine;
            op_on_machine[i].emplace_back(machine);
            std::cin >> P[i][j][machine].l >> P[i][j][machine].m >> P[i][j][machine].r;
            //std::cout << P[i][j][machine].l << " " << P[i][j][machine].m << " " << P[i][j][machine].r << '\n';
        }
    }
    for (int i = 1; i <= n; ++i)
    {
        std::cin >> D[i].m >> D[i].r;
        D[i].l = 0;
    }
    
    // DE
    while (itn--)
    {
        init();
        mutation();
        crossover();
        selection();
        //improved_selection();
        //std::cout << best_fitness << '\n';
    }
    std::cout << best_fitness << '\n';
    print_sol(best_sol);

    return 0;
}
