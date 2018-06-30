#ifndef __ITERATIVE_METHOD_HPP__
#define __ITERATIVE_METHOD_HPP__

#include <cmath>

#include "./matrix/matrix.hpp"

#include "lup.hpp"

typedef unsigned order;
static constexpr const order order_infinity = std::numeric_limits<order>::max();

typedef unsigned iteration_count;

template <typename T, order P>
Matrix<T> successive_iteration(const Matrix<T> &, const Matrix<T> &,
    typename Matrix<T>::floating_point, iteration_count &, iteration_count &);

template <typename T, order P>
Matrix<T> gauss_seidel(const Matrix<T> &a, const Matrix<T> &b,
    typename Matrix<T>::floating_point, iteration_count &, iteration_count &);

template <typename T, order P>
class Mu
{
public:
    typename Matrix<T>::floating_point operator ()(const Matrix<T> &) const;
};

template <typename T>
class Mu<T, 1u>
{
public:
    typename Matrix<T>::floating_point operator ()(const Matrix<T> &) const;
};

template <typename T>
class Mu<T, order_infinity>
{
public:
    typename Matrix<T>::floating_point operator ()(const Matrix<T> &) const;
};

template <typename T>
void arg_vilid(typename Matrix<T>::size_type, typename Matrix<T>::size_type,
    typename Matrix<T>::size_type, typename Matrix<T>::size_type,
    typename Matrix<T>::floating_point);

template <typename T>
void equivalent_form(Matrix<T> &, Matrix<T> &);

template <typename T, order P>
Matrix<T> gauss_seidel_arg(const Matrix<T> &, const Matrix<T> &, Matrix<T> &, Matrix<T> &,
    const Mu<T, P> &,
    typename Matrix<T>::floating_point &, typename Matrix<T>::floating_point &,
    typename Matrix<T>::floating_point &);

template <typename T>
iteration_count theoretical_k(typename Matrix<T>::floating_point, typename Matrix<T>::floating_point,
    typename Matrix<T>::floating_point);

template <typename T, order P>
typename Matrix<T>::floating_point epsilon_k(const Mu<T, P> &, typename Matrix<T>::floating_point,
    const Matrix<T> &, const Matrix<T> &);

template <typename T, order P>
Matrix<T> successive_iteration(const Matrix<T> &a, const Matrix<T> &b,
    const typename Matrix<T>::floating_point epsilon,
    iteration_count &th_k, iteration_count &r_k)
{
    arg_vilid<T>(a.size1(), a.size2(), b.size1(), b.size2(), epsilon);
    Matrix<T> alpha = a, beta = b;
    equivalent_form(alpha, beta);

    const Mu<T, P> mu;
    const typename Matrix<T>::floating_point mu_alpha = mu(alpha), mu_beta = mu(beta);
    Matrix<T> x_k = beta, x_k_minus_1;
    typename Matrix<T>::floating_point epsilon_k_th;
    r_k = 0u;
    do
    {
        x_k_minus_1 = std::move(x_k);
        x_k = beta + alpha * x_k_minus_1;
        epsilon_k_th = epsilon_k(mu, mu_alpha, x_k, x_k_minus_1);
        ++r_k;
    } while (std::isgreater(epsilon_k_th, epsilon));

    th_k = theoretical_k<T>(mu_alpha, mu_beta, epsilon);
    return x_k;
}

template <typename T, order P>
Matrix<T> gauss_seidel(const Matrix<T> &a, const Matrix<T> &b,
    const typename Matrix<T>::floating_point epsilon,
    iteration_count &th_k, iteration_count &r_k)
{
    arg_vilid<T>(a.size1(), a.size2(), b.size1(), b.size2(), epsilon);

    Matrix<T> first, second, x_k, x_k_minus_1;
    const Mu<T, P> mu;
    typename Matrix<T>::floating_point mu_alpha, mu_beta, mu_c;
    x_k = gauss_seidel_arg(a, b, first, second, mu, mu_alpha, mu_beta, mu_c);

    typename Matrix<T>::floating_point epsilon_k_th;
    r_k = 0u;
    do
    {
        x_k_minus_1 = std::move(x_k);
        x_k = second + first * x_k_minus_1;
        epsilon_k_th = epsilon_k(mu, mu_alpha, x_k, x_k_minus_1);
        ++r_k;
    } while (std::isgreater(epsilon_k_th, epsilon));

    th_k = theoretical_k<T>(mu_alpha, mu_beta, epsilon);
    return x_k;
}

template <typename T, order P>
typename Matrix<T>::floating_point Mu<T, P>::operator ()(const Matrix<T> &a) const
{
    if (!a.size1() || !a.size2())
    {
        throw std::domain_error("Can't calculate");
    }

    const typename Matrix<T>::floating_point p = P;
    typename Matrix<T>::floating_point result = 0.0;
    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
        {
            result += std::pow(std::fabs(a(i, j)), p);
        }
    }

    return std::pow(result, 1.0 / p);
}


template <typename T>
typename Matrix<T>::floating_point Mu<T, 1u>::operator ()(const Matrix<T> &a) const
{
    if (!a.size1() || !a.size2())
    {
        throw std::domain_error("Can't calculate");
    }

    typename Matrix<T>::floating_point result = -1.0;
    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
    {
        typename Matrix<T>::floating_point column_mu = 0.0;
        for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
        {
            column_mu += std::fabs(a(i, j));
        }
        if (std::isless(result, column_mu))
        {
            result = column_mu;
        }
    }

    return result;
}

template <typename T>
typename Matrix<T>::floating_point Mu<T, order_infinity>::operator ()(const Matrix<T> &a) const
{
    if (!a.size1() || !a.size2())
    {
        throw std::domain_error("Can't calculate");
    }

    typename Matrix<T>::floating_point result = -1.0;
    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        typename Matrix<T>::floating_point row_mu = 0.0;
        for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
        {
            row_mu += std::fabs(a(i, j));
        }
        if (std::isless(result, row_mu))
        {
            result = row_mu;
        }
    }

    return result;
}

template <typename T>
void arg_vilid(const typename Matrix<T>::size_type a1, const typename Matrix<T>::size_type a2,
    const typename Matrix<T>::size_type b1, const typename Matrix<T>::size_type b2,
    const typename Matrix<T>::floating_point epsilon)
{
    if (a1 != a2 || !a2 || a1 != b1 || b2 != 1u ||
        std::isless(epsilon, Matrix<T>::value_epsilon))
    {
        throw std::domain_error("Can't be solved");
    }
}

template <typename T>
void equivalent_form(Matrix<T> &alpha, Matrix<T> &beta)
{
    const Matrix<T> p = LUP<T>(alpha).get_p();
    alpha = p * alpha, beta = p * beta;
    const typename Matrix<T>::size_type m = alpha.size1(), n = alpha.size2();
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        const typename Matrix<T>::floating_point divider = alpha(i, i);
        for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
        {
            if (j == i)
            {
                alpha(i, j) = 0.0;
            }
            else
            {
                alpha(i, j) /= -divider;
            }
        }
        beta(i) /= divider;
    }
}

template <typename T, order P>
Matrix<T> gauss_seidel_arg(const Matrix<T> &a, const Matrix<T> &b, Matrix<T> &first, Matrix<T> &second,
    const Mu<T, P> &mu,
    typename Matrix<T>::floating_point &mu_alpha,
    typename Matrix<T>::floating_point &mu_beta,
    typename Matrix<T>::floating_point &mu_c)
{
    Matrix<T> alpha = a, beta = b;
    equivalent_form(alpha, beta);
    mu_alpha = mu(alpha);
    mu_beta = mu(beta);

    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    Matrix<T> b_matrix(m, n), c_matrix(m, n);
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
        {
            if (i < j)
            {
                b_matrix(i, j) = alpha(i, j);
            }
            else
            {
                c_matrix(i, j) = alpha(i, j);
            }
        }
    }
    mu_c = mu(c_matrix);

    Matrix<T> invertible = LUP<T>(Matrix<T>().identity(m) - b_matrix).invertible();
    first = invertible * c_matrix;
    second = invertible * beta;

    return beta;
}

template <typename T>
iteration_count theoretical_k(const typename Matrix<T>::floating_point mu_alpha,
    const typename Matrix<T>::floating_point mu_beta,
    const typename Matrix<T>::floating_point epsilon)
{
    return std::ceil((std::log(epsilon) + std::log(1.0 - mu_alpha) - std::log(mu_beta)) / std::log(mu_alpha) - 1.0);
}

template <typename T, order P>
typename Matrix<T>::floating_point epsilon_k(const Mu<T, P> &mu, typename Matrix<T>::floating_point mu_alpha,
    const Matrix<T> &x_k, const Matrix<T> &x_k_minus_1)
{
    typename Matrix<T>::floating_point result = mu(x_k - x_k_minus_1);
    if (mu_alpha < 1.0)
    {
        result *= mu_alpha / (1.0 - mu_alpha);
    }

    return result;
}

#endif
