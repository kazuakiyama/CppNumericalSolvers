/**
 * Copyright (c) 2014-2015 Patrick Wieschollek
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "LbfgsbSolver.hpp"
#include <iostream>
#include <list>
#include "stopwatch.hpp"

namespace pwie
{

template <typename _Scalar>
std::vector<int>
sort_indexes(const std::vector< std::pair<int, _Scalar> > & v)
{
    std::vector<int> idx(v.size());
    for(size_t i = 0; i != idx.size(); ++i)
        idx[i] = v[i].first;
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2)
    {
        return v[i1].second < v[i2].second;
    });
    return idx;
}


template <typename Func>
LbfgsbSolver<Func>::LbfgsbSolver(const Func & func) : ISolver<Func>(func)
{
    // TODO Auto-generated constructor stub
}

template <typename Func>
void
LbfgsbSolver<Func>::GetGeneralizedCauchyPoint(const InputType & x, JacobianType & g,
                                              InputType & x_cauchy, VectorX & c)
{
    // PAGE 8
    // Algorithm CP: Computation of the generalized Cauchy point
    // Given x,l,u,g, and B = \theta I-WMW

    // {all t_i} = { (idx,value), ... }
    // TODO: use "std::set" ?
    std::vector<std::pair<int, Scalar> > SetOfT;
    // the feasible set is implicitly given by "SetOfT - {t_i==0}"
    JacobianType d = JacobianType::Zero(_DIM, 1);

    // n operations
    for(int j = 0; j < _DIM; j++)
    {
        if(g(j) == 0)
        {
            SetOfT.push_back(std::make_pair(j, std::numeric_limits<Scalar>::max()));
        }
        else
        {
            Scalar tmp = 0;
            if(g(j) < 0)
            {
                tmp = (x(j) - _ub(j)) / g(j);
            }
            else
            {
                tmp = (x(j) - _lb(j)) / g(j);
            }
            d(j) = -g(j);
            SetOfT.push_back(std::make_pair(j, tmp));
        }

    }
    Debug("GGCP" << d.transpose());

    // paper: using heapsort
    // sortedindices [1,0,2] means the minimal element is on the 1th entry
    std::vector<int> SortedIndices = sort_indexes(SetOfT);

    x_cauchy = x;
    // Initialize
    // p :=     W^T*p
    VectorX p = (W.transpose() * d);                     // (2mn operations)
    // c :=     0
    c = MatrixX::Zero(M.rows(), 1);
    // f' :=    g^T*d = -d^Td
    Scalar f_prime = -d.dot(d);                         // (n operations)
    // f'' :=   \theta*d^T*d-d^T*W*M*W^T*d = -\theta*f' - p^T*M*p
    Scalar f_doubleprime = (Scalar)(-1.0 * theta) * f_prime - p.dot(M * p); // (O(m^2) operations)
    // \delta t_min :=  -f'/f''
    Scalar dt_min = -f_prime / f_doubleprime;
    // t_old :=     0
    Scalar t_old = 0;
    // b :=     argmin {t_i , t_i >0}
    int i = 0;
    for(int j = 0; j < _DIM; j++)
    {
        i = j;
        if(SetOfT[SortedIndices[j]].second != 0)
            break;
    }
    int b = SortedIndices[i];
    // see below
    // t                    :=  min{t_i : i in F}
    Scalar t = SetOfT[b].second;
    // \delta t             :=  t - 0
    Scalar dt = t - t_old;

    // examination of subsequent segments
    while((dt_min >= dt) && (i < _DIM))
    {
        if(d(b) > 0)
            x_cauchy(b) = _ub(b);
        else if(d(b) < 0)
            x_cauchy(b) = _lb(b);

        // z_b = x_p^{cp} - x_b
        Scalar zb = x_cauchy(b) - x(b);
        // c   :=  c +\delta t*p
        c += dt * p;
        // cache
        VectorX wbt = W.row(b);
        f_prime += dt * f_doubleprime + (Scalar) g(b) * g(b)
                   + (Scalar) theta * g(b) * zb
                   - (Scalar) g(b) * wbt.transpose() * (M * c);
        f_doubleprime += (Scalar) - 1.0 * theta * g(b) * g(b)
                         - (Scalar) 2.0 * (g(b) * (wbt.dot(M * p)))
                         - (Scalar) g(b) * g(b) * wbt.transpose() * (M * wbt);
        p += g(b) * wbt.transpose();
        d(b) = 0;
        dt_min = -f_prime / f_doubleprime;
        t_old = t;
        ++i;
        if(i < _DIM)
        {
            b = SortedIndices[i];
            t = SetOfT[b].second;
            dt = t - t_old;
        }

    }

    dt_min = MAX(dt_min, 0);
    t_old += dt_min;

    Debug(SortedIndices[0] << " " << SortedIndices[1]);

#pragma omp parallel for
    for(int ii = i; ii < x_cauchy.rows(); ii++)
    {
        x_cauchy(SortedIndices[ii]) = x(SortedIndices[ii])
                                      + t_old * d(SortedIndices[ii]);
    }
    Debug(x_cauchy.transpose());

    c += dt_min * p;
    Debug(c.transpose());
}
template <typename Func>
typename LbfgsbSolver<Func>::Scalar
LbfgsbSolver<Func>::FindAlpha(const InputType & x_cp, const VectorX & du, const std::vector<int> & FreeVariables)
{
    /* this returns
     * a* = max {a : a <= 1 and  l_i-xc_i <= a*d_i <= u_i-xc_i}
     */
    Scalar alphastar = 1;
    const unsigned int n = FreeVariables.size();
    for(unsigned int i = 0; i < n; i++)
    {
        if(du(i) > 0)
        {
            alphastar = MIN(alphastar,
                            (_ub(FreeVariables[i]) - x_cp(FreeVariables[i]))
                            / du(i));
        }
        else
        {
            alphastar = MIN(alphastar,
                            (_lb(FreeVariables[i]) - x_cp(FreeVariables[i]))
                            / du(i));
        }
    }
    return alphastar;
}

template <typename Func>
void
LbfgsbSolver<Func>::SubspaceMinimization(InputType & x_cauchy, InputType & x, const VectorX & c,
                                         JacobianType & g, InputType & SubspaceMin)
{
    // cached value: ThetaInverse=1/theta;
    Scalar theta_inverse = 1 / theta;

    // size of "t"
    std::vector<int> FreeVariablesIndex;
    Debug(x_cauchy.transpose());

    //std::cout << "free vars " << FreeVariables.rows() << std::endl;
    for(int i = 0; i < x_cauchy.rows(); i++)
    {
        Debug(x_cauchy(i) << " " << _ub(i) << " " << _lb(i));
        if((x_cauchy(i) != _ub(i)) && (x_cauchy(i) != _lb(i)))
        {
            FreeVariablesIndex.push_back(i);
        }
    }
    const int FreeVarCount = FreeVariablesIndex.size();

    MatrixX WZ = MatrixX::Zero(W.cols(), FreeVarCount);

    for(int i = 0; i < FreeVarCount; i++)
        WZ.col(i) = W.row(FreeVariablesIndex[i]);

    Debug(WZ);

    // r=(g+theta*(x_cauchy-x)-W*(M*c));
    Debug(g);
    Debug(x_cauchy);
    Debug(x);
    InputType rr = (g + theta * (x_cauchy - x) - W * (M * c));
    // r=r(FreeVariables);
    VectorX r = VectorX::Zero(FreeVarCount, 1);
    for(int i = 0; i < FreeVarCount; i++)
        r.row(i) = rr.row(FreeVariablesIndex[i]);

    Debug(r.transpose());

    // STEP 2: "v = w^T*Z*r" and STEP 3: "v = M*v"
    VectorX v = M * (WZ * r);
    // STEP 4: N = 1/theta*W^T*Z*(W^T*Z)^T
    MatrixX N = theta_inverse * WZ * WZ.transpose();
    // N = I - MN
    N = MatrixX::Identity(N.rows(), N.rows()) - M * N;
    // STEP: 5
    // v = N^{-1}*v
    v = N.lu().solve(v);
    // STEP: 6
    // HERE IS A MISTAKE IN THE ORIGINAL PAPER!
    VectorX du = -theta_inverse * r
                  - theta_inverse * theta_inverse * WZ.transpose() * v;
    Debug(du.transpose());
    // STEP: 7
    Scalar alpha_star = FindAlpha(x_cauchy, du, FreeVariablesIndex);

    // STEP: 8
    VectorX dStar = alpha_star * du;

    SubspaceMin = x_cauchy;
    for(int i = 0; i < FreeVarCount; i++)
    {
        SubspaceMin(FreeVariablesIndex[i]) = SubspaceMin(
                FreeVariablesIndex[i]) + dStar(i);
    }
}

template <typename Func>
void
LbfgsbSolver<Func>::internalSolve(InputType & x0)
{
    _DIM = x0.rows();

    _lb = _functor.getLowerBound(_DIM);
    _ub = _functor.getUpperBound(_DIM);
    theta = 1.0;

    W = MatrixX::Zero(_DIM, 0);
    M = MatrixX::Zero(0, 0);

    Assert(x0.rows() == _lb.rows(), "lower bound size incorrect");
    Assert(x0.rows() == _ub.rows(), "upper bound size incorrect");

    Debug(x0.transpose());
    Debug(_lb.transpose());
    Debug(_ub.transpose());

    Assert((x0.array() >= _lb.array()).all(),
           "seed is not feasible (violates lower bound)");
    Assert((x0.array() <= _ub.array()).all(),
           "seed is not feasible (violates upper bound)");

    MatrixX yHistory = MatrixX::Zero(_DIM, 0);
    MatrixX sHistory = MatrixX::Zero(_DIM, 0);

    InputType x = x0;
    JacobianType g(_DIM);
    size_t k = 0;

    Scalar f = _functor.f(x);
    _functor.gradient(x, g);

    Debug(f);
    Debug(g.transpose());

    auto noConvergence =
        [&](InputType & x, JacobianType & g)->bool
    {
      auto d = (x - g).cwiseMax(_lb).cwiseMin(_ub) - x;
      return d.template lpNorm<Eigen::Infinity>() >= settings.tol;
    };

    int no_improve_count = 0;
    Stopwatch<> stopwatch;
    while(noConvergence(x, g) && (k < settings.maxIter))
    {
        Scalar f_old = f;
        InputType x_old = x;
        JacobianType g_old = g;

        // STEP 2: compute the cauchy point by algorithm CP
        InputType CauchyPoint = MatrixX::Zero(_DIM, 1);
        VectorX c = MatrixX::Zero(_DIM, 1);
        GetGeneralizedCauchyPoint(x, g, CauchyPoint, c);
        // STEP 3: compute a search direction d_k by the primal method
        InputType SubspaceMin;
        SubspaceMinimization(CauchyPoint, x, c, g, SubspaceMin);

        Scalar step = 2;

        // STEP 4: perform linesearch and STEP 5: compute gradient
        int searchstatus = LineSearch(x, SubspaceMin - x, f, g, step);
        if (searchstatus < 0) {
          return;
        }

        // prepare for next iteration
        JacobianType newY = g - g_old;
        InputType newS = x - x_old;

        // STEP 6:
        Scalar test = newS.dot(newY);
        test = (test < 0) ? -1.0 * test : test;

        if(test > (std::numeric_limits<Scalar>::min() * 1e-6 * newY.squaredNorm()))
        {
            if(k < settings.m)
            {
                yHistory.conservativeResize(_DIM, k + 1);
                sHistory.conservativeResize(_DIM, k + 1);
            }
            else
            {
                yHistory.leftCols(settings.m - 1) = 
                  yHistory.rightCols(settings.m - 1).eval();
                sHistory.leftCols(settings.m - 1) = 
                  sHistory.rightCols(settings.m - 1).eval();
            }
            yHistory.rightCols(1) = newY;
            sHistory.rightCols(1) = newS;

            // STEP 7:
            theta = (Scalar)(newY.transpose() * newY)
                    / (newY.transpose() * newS);

            W = MatrixX::Zero(yHistory.rows(),
                              yHistory.cols() + sHistory.cols());

            W << yHistory, (theta * sHistory);

            MatrixX A = sHistory.transpose() * yHistory;
            MatrixX L = A.template triangularView<Eigen::StrictlyLower>();
            MatrixX MM(A.rows() + L.rows(), A.rows() + L.cols());
            MatrixX D = -1 * A.diagonal().asDiagonal();
            MM << D, L.transpose(), L, ((sHistory.transpose() * sHistory)
                                        * theta);

            M = MM.inverse();
        }

        Scalar ttt;
        ttt = f_old - f;
        using std::abs;
        Debug("--> " << abs(ttt));
        if (abs(ttt) < 1e-13)
        {
          // successive function values too similar
          if (no_improve_count++ > 10) {
            std::cout << "no improvement, norm=" << abs(ttt) << "\n";
            break;
          }
        }
        else {
          no_improve_count = 0;
        }

        k++;
        settings.numIters = k;

        InputType dir = SubspaceMin - x;
        ISolver<Func>::checkConverged(k, x, step, dir, g);

        if (settings.verbosity > 0) {
          stopwatch.stop();
          std::cout << "iteration " << k << " " << f
                    << " dt=" << stopwatch.elapsed()/1e3
                    << " df=" << ttt
                    << " step=" << step << std::endl;
          stopwatch.start();
        }
    }

    x0 = x;
}
}

/* namespace pwie */
