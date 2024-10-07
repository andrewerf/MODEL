//
// Created by Andrey Aralov on 10/7/24.
//
#pragma once
#include <iostream>

template <typename T, Index n, Index m>
void print( const Mat<T, n, m>& mat )
{
    for ( Index i = 0; i < n; ++i )
    {
        for ( Index j = 0; j < m; ++j )
        {
            std::cout << mat( i, j ) << ' ';
        }
        std::cout << '\n';
    }
}