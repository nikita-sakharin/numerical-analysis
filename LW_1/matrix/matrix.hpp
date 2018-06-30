#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <type_traits>

template <typename T>
class Matrix
{
public:
    typedef typename std::enable_if<std::is_floating_point<T>::value, T>::type floating_point;
    typedef floating_point value_type;

private:
    typedef std::vector<value_type> row_type;

public:
    typedef typename row_type::reference reference;
    typedef typename row_type::const_reference const_reference;
    typedef typename row_type::size_type size_type;
    typedef typename row_type::difference_type difference_type;

    static constexpr const floating_point value_epsilon = std::numeric_limits<floating_point>::epsilon();
    static constexpr const floating_point value_round_error = std::numeric_limits<floating_point>::round_error();
    static constexpr const size_type size_min = std::numeric_limits<size_type>::min();
    static constexpr const size_type size_max = std::numeric_limits<size_type>::max();

private:
    std::vector<row_type> matrix;

public:
    Matrix() noexcept = default;
    Matrix(size_type, size_type);
    Matrix(size_type, size_type, const value_type &);
    Matrix(const Matrix &) = default;
    Matrix(Matrix &&) noexcept = default;
    Matrix &operator =(const Matrix &) = default;
    Matrix &operator =(Matrix &&) = default;
    ~Matrix() = default;

    size_type size1() const noexcept;
    size_type size2() const noexcept;
    void resize(size_type, size_type);
    reference operator ()(size_type, size_type);
    const_reference operator ()(size_type, size_type) const;
    reference operator ()(size_type);
    const_reference operator ()(size_type) const;

    Matrix &zero(size_type, size_type);
    Matrix &identity(size_type);

    template <typename S>
    friend Matrix<S> transpose(const Matrix<S> &);
    Matrix &transpose();

    void row_switching(size_type, size_type);
    void row_multiplication(value_type, size_type);
    void row_addition(size_type, value_type, size_type);

    Matrix operator +() const;
    Matrix operator -() const;

    template <typename S>
    friend Matrix<S> operator *(const Matrix<S> &, const typename Matrix<S>::value_type &);
    template <typename S>
    friend Matrix<S> operator *(const typename Matrix<S>::value_type &, const Matrix<S> &);
    Matrix &operator *=(const value_type &) noexcept;
    template <typename S>
    friend Matrix<S> operator /(const Matrix<S> &, const typename Matrix<S>::value_type &);
    Matrix &operator /=(const value_type &) noexcept;

    template <typename S>
    friend Matrix<S> operator +(const Matrix<S> &, const Matrix<S> &);
    Matrix &operator +=(const Matrix &);

    template <typename S>
    friend Matrix<S> operator -(const Matrix<S> &, const Matrix<S> &);
    Matrix &operator -=(const Matrix &);

    template <typename S>
    friend Matrix<S> operator *(const Matrix<S> &, const Matrix<S> &);
    Matrix &operator *=(const Matrix &);

    template <typename S>
    friend std::istream &operator >>(std::istream &, Matrix<S> &);
    template <typename S>
    friend std::ostream &operator <<(std::ostream &, const Matrix<S> &);
};

template <typename T>
Matrix<T>::Matrix(const Matrix<T>::size_type row,
    const Matrix<T>::size_type column) : matrix(row, row_type(column))
{
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>::size_type row,
    const Matrix<T>::size_type column,
    const Matrix<T>::value_type &cell) : matrix(row, row_type(column, cell))
{
}

template <typename T>
typename Matrix<T>::size_type Matrix<T>::size1() const noexcept
{
    return matrix.size();
}

template <typename T>
typename Matrix<T>::size_type Matrix<T>::size2() const noexcept
{
    if (matrix.size())
    {
        return matrix.front().size();
    }
    else
    {
        return 0u;
    }
}

template <typename T>
void Matrix<T>::resize(const Matrix<T>::size_type row_cnt, const Matrix<T>::size_type column_cnt)
{
    matrix.resize(row_cnt);
    for (row_type &row : matrix)
    {
        row.resize(column_cnt);
    }
}

template <typename T>
typename Matrix<T>::reference Matrix<T>::operator ()(const Matrix<T>::size_type i,
    const Matrix<T>::size_type j)
{
    return matrix.at(i).at(j);
}

template <typename T>
typename Matrix<T>::const_reference Matrix<T>::operator ()(const Matrix<T>::size_type i,
    const Matrix<T>::size_type j) const
{
    return matrix.at(i).at(j);
}

template <typename T>
typename Matrix<T>::reference Matrix<T>::operator ()(const Matrix<T>::size_type i)
{
    if (matrix.size() == 1u)
    {
        return matrix[0u].at(i);
    }
    else if (matrix.at(0u).size() == 1u)
    {
        return matrix.at(i)[0u];
    }
    else
    {
        throw std::domain_error("Is't vector");
    }
}

template <typename T>
typename Matrix<T>::const_reference Matrix<T>::operator ()(const Matrix<T>::size_type i) const
{
    if (matrix.size() == 1u)
    {
        return matrix[0u].at(i);
    }
    else if (matrix.at(0u).size() == 1u)
    {
        return matrix.at(i)[0u];
    }
    else
    {
        throw std::domain_error("Is't vector");
    }
}

template <typename T>
Matrix<T> &Matrix<T>::zero(const Matrix<T>::size_type row_cnt,
    const Matrix<T>::size_type column_cnt)
{
    matrix.assign(row_cnt, row_type(column_cnt));
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::identity(const Matrix<T>::size_type order)
{
    matrix.assign(order, row_type(order));
    for (size_type i = 0u; i < order; ++i)
    {
        matrix[i][i] = 1.0;
    }

    return *this;
}

template <typename T>
Matrix<T> transpose(const Matrix<T> &rhs)
{
    const typename Matrix<T>::size_type m = rhs.matrix.size(), n = rhs.matrix.front().size();
    Matrix<T> result(n, m);
    for (typename Matrix<T>::size_type i = 0u; i < n; ++i)
    {
        for (typename Matrix<T>::size_type j = 0u; j < m; ++j)
        {
            result.matrix[i][j] = rhs.matrix[j][i];
        }
    }

    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::transpose()
{
    const size_type m = matrix.size(), n = matrix.front().size();
    if (m == n)
    {
        for (size_type i = 1u; i < m; ++i)
        {
            for (size_type j = 0u; j < i; ++j)
            {
                std::swap(matrix[i][j], matrix[j][i]);
            }
        }
    }
    else
    {
        *this = transpose(*this);
    }

    return *this;
}

template <typename T>
void Matrix<T>::row_switching(const Matrix<T>::size_type i,
    const Matrix<T>::size_type j)
{
    std::swap(matrix[i], matrix[j]);
}

template <typename T>
void Matrix<T>::row_multiplication(const Matrix<T>::value_type alpha,
    const Matrix<T>::size_type i)
{
    const size_type size = matrix.front().size();
    for (size_type k = 0u; k < size; ++k)
    {
        matrix[i][k] *= alpha;
    }
}

template <typename T>
void Matrix<T>::row_addition(const Matrix<T>::size_type i,
    const Matrix<T>::value_type alpha,
    const Matrix<T>::size_type j)
{
    const size_type size = matrix.front().size();
    for (size_type k = 0u; k < size; ++k)
    {
        matrix[i][k] = std::fma(alpha, matrix[j][k], matrix[i][k]);
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator +() const
{
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator -() const
{
    Matrix<T> result = *this;
    const size_type m = matrix.size(), n = matrix.front().size();
    for (size_type i = 0u; i < m; ++i)
    {
        for (size_type j = 0u; j < n; ++j)
        {
            result.matrix[i][j] = -matrix[i][j];
        }
    }

    return result;
}

template <typename T>
Matrix<T> operator *(const Matrix<T> &lhs, const typename Matrix<T>::value_type &value)
{
    Matrix<T> result = lhs;
    result *= value;
    return result;
}

template <typename T>
Matrix<T> operator *(const typename Matrix<T>::value_type &value, const Matrix<T> &rhs)
{
    return rhs * value;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *=(const Matrix<T>::value_type &value) noexcept
{
    const size_type m = matrix.size(), n = matrix.front().size();
    for (size_type i = 0u; i < m; ++i)
    {
        for (size_type j = 0u; j < n; ++j)
        {
            matrix[i][j] *= value;
        }
    }

    return *this;
}

template <typename T>
Matrix<T> operator /(const Matrix<T> &lhs, const typename Matrix<T>::value_type &value)
{
    Matrix<T> result = lhs;
    result /= value;
    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::operator /=(const Matrix<T>::value_type &value) noexcept
{
    if (std::isless(std::fabs(value), value_epsilon))
    {
        throw std::overflow_error("Division by zero");
    }
    const size_type m = matrix.size(), n = matrix.front().size();
    for (size_type i = 0u; i < m; ++i)
    {
        for (size_type j = 0u; j < n; ++j)
        {
            matrix[i][j] /= value;
        }
    }

    return *this;
}

template <typename T>
Matrix<T> operator +(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
    Matrix<T> result = lhs;
    result += rhs;
    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::operator +=(const Matrix<T> &rhs)
{
    const size_type m = matrix.size(), n = matrix.front().size();
    if (m != rhs.matrix.size() || n != rhs.matrix.front().size())
    {
        throw std::domain_error("Matrices can't be summed");
    }

    for (size_type i = 0u; i < m; ++i)
    {
        for (size_type j = 0u; j < n; ++j)
        {
            matrix[i][j] += rhs.matrix[i][j];
        }
    }

    return *this;
}

template <typename T>
Matrix<T> operator -(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
    Matrix<T> result = lhs;
    result -= rhs;
    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::operator -=(const Matrix<T> &rhs)
{
    const size_type m = matrix.size(), n = matrix.front().size();
    if (m != rhs.matrix.size() || n != rhs.matrix.front().size())
    {
        throw std::domain_error("Matrices can't be subtracted");
    }

    for (size_type i = 0u; i < m; ++i)
    {
        for (size_type j = 0u; j < n; ++j)
        {
            matrix[i][j] -= rhs.matrix[i][j];
        }
    }

    return *this;
}

template <typename T>
Matrix<T> operator *(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
    const typename Matrix<T>::size_type m = lhs.matrix.size(), l = lhs.matrix.front().size(), n = rhs.matrix.front().size();
    if (l != rhs.matrix.size())
    {
        throw std::domain_error("Matrices can't be multiplied");
    }

    Matrix<T> result;
    result.zero(m, n);
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        for (typename Matrix<T>::size_type j = 0u; j < n; ++j)
        {
            for (typename Matrix<T>::size_type k = 0u; k < l; ++k)
            {
                result.matrix[i][j] = std::fma(lhs.matrix[i][k], rhs.matrix[k][j], result.matrix[i][j]);
            }
        }
    }

    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *=(const Matrix<T> &rhs)
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
std::istream &operator >>(std::istream &is, Matrix<T> &m)
{
    typename Matrix<T>::size_type row_cnt, column_cnt;
    is >> row_cnt >> column_cnt;
    m.matrix.assign(row_cnt, typename Matrix<T>::row_type(column_cnt));
    for (typename Matrix<T>::row_type &row : m.matrix)
    {
        for (typename Matrix<T>::value_type &value : row)
        {
            is >> value;
        }
    }

    return is;
}

template <typename T>
std::ostream &operator <<(std::ostream &os, const Matrix<T> &m)
{
    os << m.matrix.size() << ' ' << m.matrix.front().size() << '\n';
    for (const typename Matrix<T>::row_type &row : m.matrix)
    {
        for (const typename Matrix<T>::value_type &value : row)
        {
            os << value << ' ';
        }
        os << '\n';
    }

    return os;
}

#endif
