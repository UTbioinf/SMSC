#ifndef __MULTI_ARRAY_H
#define __MULTI_ARRAY_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace loon
{

/* ====================== class Array2D ==========================
================================================================ */
template<class T>
class Array2D: public std::vector<T>
{
private:
    size_t column;
    size_t row;

public:
    Array2D(size_t n = 0, size_t c=0, size_t r = 0);
    void set_column(size_t c);
    void resize(size_t row);
    void resize(size_t row, size_t col); // this will also change the column of the 2-D array
    void resize_column(size_t c);
    void add_row();
    T& at(size_t i, size_t j);
    const T& at(size_t i, size_t j) const;
    void clear();
    size_t get_row() const;
};

template<class T>
Array2D<T>::Array2D(size_t n /*=0*/, size_t c /*=0*/, size_t r /*=0*/): 
        column(c), row(r), 
        std::vector<T>::vector(n)
{
}

template<class T>
void Array2D<T>::set_column(size_t c)
{
    column = c;
}

template<class T>
void Array2D<T>::resize(size_t row)
{
    this->row = row;
    if(row * column > this->size())
        std::vector<T>::resize( std::max(row * column, this->size() << 1) );
}

template<class T>
void Array2D<T>::resize(size_t row, size_t col)
{
    column = col;
    this->row = row;
    if(row * col > this->size())
        std::vector<T>::resize( std::max(row * col, this->size() << 1) );
}

template<class T>
void Array2D<T>::resize_column(size_t c)
{
    row = 1;
    column = c;
    if(column > this->size())
        std::vector<T>::resize( std::max(column, this->size() << 1) );
}

template<class T>
void Array2D<T>::add_row()
{
    ++row;
    if(row * column > this->size())
        std::vector<T>::resize( std::max(row * column, this->size() << 1) );
}

template<class T>
T& Array2D<T>::at(size_t i, size_t j)
{
    if(row <= i || column <= j)
        throw std::out_of_range("Array2D_range_check");
    return std::vector<T>::at( i * column + j);
}

template<class T>
const T& Array2D<T>::at(size_t i, size_t j) const
{
    if(row <= i || column <= j)
        throw std::out_of_range("Array2D_range_check");
    return std::vector<T>::at( i * column + j);
}

template<class T>
void Array2D<T>::clear()
{
    column = row = 0;
    std::vector<T>::clear();
}

template<class T>
size_t Array2D<T>::get_row() const
{
    return row;
}


/* ====================== class Square2D ==========================
================================================================ */
template<class T>
class Square2D: public std::vector<T>
{
private:
    size_t n;
    T default_value;
public:
    Square2D();
    Square2D(size_t reserve, const T& def_val = T());
    void set_default(const T& def_val);
    void reset();
    T& at(size_t i, size_t j);
    const T& at(size_t i, size_t j) const;
    void setn(size_t nn);
    size_t getn() const;
};

template<class T>
Square2D<T>::Square2D(): n(0), default_value(0)
{
}

template<class T>
Square2D<T>::Square2D(size_t reserve, const T& def_val/* = T() */):
        n(0), default_value(def_val), std::vector<T>::vector(n, def_val)
{
}

template<class T>
void Square2D<T>::set_default(const T& def_val)
{
    default_value = def_val;
}

template<class T>
void Square2D<T>::reset()
{
    n = 0;
}

template<class T>
T& Square2D<T>::at(size_t i, size_t j)
{
    if(i < j)
    {
        if(n <= j)
        {
            size_t old_cap = n * n;
            n = j + 1;
            size_t new_cap = n * n;
            if(std::vector<T>::size() < new_cap)
                std::vector<T>::resize( std::max(std::vector<T>::size() << 2, new_cap));
            for(size_t ii = old_cap; ii < new_cap; ++ii)
                std::vector<T>::at(ii) = default_value;
        }
        return std::vector<T>::at( j*j + i );
    }
    else
    {
        if(n <= i)
        {
            size_t old_cap = n * n;
            n = i + 1;
            size_t new_cap = n * n;
            if(std::vector<T>::size() < new_cap)
                std::vector<T>::resize( std::max(std::vector<T>::size() << 2, new_cap) );
            for(size_t ii = old_cap; ii < new_cap; ++ii)
                std::vector<T>::at(ii) = default_value;
        }
        return std::vector<T>::at( i*i + i + j );
    }
}

template<class T>
const T& Square2D<T>::at(size_t i, size_t j) const
{
    if(n <= std::max(i, j))
        throw std::out_of_range("Square2D_range_check");
    if(i < j)
        return std::vector<T>::at( j*j + i);
    else
        return std::vector<T>::at( i*i + i + j );
}

template<class T>
void Square2D<T>::setn(size_t nn)
{
    size_t new_cap = nn * nn;
    if(n < nn)
    {
        size_t old_cap = n * n;
        if(std::vector<T>::size() < new_cap)
            std::vector<T>::resize( std::max(std::vector<T>::size() << 2, new_cap) );
    }
    std::vector<T>::assign(new_cap, default_value);
    n = nn;
}

template<class T>
size_t Square2D<T>::getn() const
{
    return n;
}

}
#endif
