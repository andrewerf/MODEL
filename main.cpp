#include <iostream>
#include <array>
#include <cassert>
#include "Mat.hpp"
#include "LU.hpp"


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





int main()
{
//    Mat<float, 2, 2> mat;
//    mat(0, 0) = 1.f;
//    mat(0, 1) = 2.f;
//    mat(1, 0) = 4.f;
//    mat(1, 1) = 5.f;

    Mat<float, 4, 4> mat;
    mat( 0, 0 ) = 1;
    mat( 0, 1 ) = 2;
    mat( 0, 2 ) = 3;
    mat( 0, 3 ) = 4;

    mat( 1, 0 ) = 1;
    mat( 1, 1 ) = 1;
    mat( 1, 2 ) = 1;
    mat( 1, 3 ) = 1;

    mat( 2, 0 ) = 4;
    mat( 2, 1 ) = 2;
    mat( 2, 2 ) = 3;
    mat( 2, 3 ) = 1;

    mat( 3, 0 ) = 2;
    mat( 3, 1 ) = 3;
    mat( 3, 2 ) = 1;
    mat( 3, 3 ) = 1;

    auto [L, U] = LU( mat );
    print( L );
    std::cout << '\n';
    print( U );
    std::cout << '\n';

    auto check = L * U;
    print( check );
    std::cout << '\n';

    return 0;
}
