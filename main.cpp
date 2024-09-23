#include <iostream>
#include <array>
#include <cassert>
#include "Mat.hpp"
#include "LU.hpp"






int main()
{
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

    auto [P, L, U, Q] = PLUQ( mat );
    print( P );
    std::cout << '\n';
    print( L );
    std::cout << '\n';
    print( U );
    std::cout << '\n';
    print( Q );
    std::cout << '\n';

    auto check = P * L * U * Q;
    print( check );
    std::cout << '\n';

//    auto [L, U] = LU( mat );
//    print( L );
//    std::cout << '\n';
//    print( U );
//    std::cout << '\n';
//
//    auto check = L * U;
//    print( check );
//    std::cout << '\n';

    return 0;
}
