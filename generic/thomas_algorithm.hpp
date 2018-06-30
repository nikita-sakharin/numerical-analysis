#ifndef __THOMAS_ALGORITHM_HPP__
#define __THOMAS_ALGORITHM_HPP__

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric;

template <typename T>
ublas::vector<T> thomas_algorithm(const ublas::vector<T> &a,
    const ublas::vector<T> &b,
    const ublas::vector<T> &c,
    const ublas::vector<T> &d)
{
    if (a.size() != d.size() || b.size() != d.size() || c.size() != d.size() ||
        d.size() < 3u)
    {
        throw std::domain_error("Can't be solved");
    }

    const typename ublas::vector<T>::size_type size = d.size();
    ublas::vector<T> p(size), q(size), x(size);

    p(0u) = -c(0u) / b(0u);
    q(0u) = d(0u) / b(0u);

    for (typename ublas::vector<T>::size_type i = 1u; i < size; ++i)
    {
        p(i) = -c(i) / (b(i) + a(i) * p(i - 1u));
        q(i) = (d(i) - a(i) * q(i - 1u)) / (b(i) + a(i) * p(i - 1u));
    }

    x(size - 1u) = q(size - 1u);

    for (typename ublas::vector<T>::difference_type i = size - 2u; i >= 0; --i)
    {
        x(i) = p(i) * x(i + 1u) + q(i);
    }

    return x;
}

#endif
