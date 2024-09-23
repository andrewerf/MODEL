//
// Created by Andrey Aralov on 9/23/24.
//
#include <iostream>
#include <array>
#include <cassert>
#include "Mat.hpp"
#include "LU.hpp"

#include <boost/multiprecision/mpfr.hpp>


void tryLU( const auto& mat )
{
    auto [L, U] = LU( mat );
    auto check = L * U;
    std::cout << "LU check:\n";
    print( check );
    std::cout << '\n';
}

void tryPLUQ( const auto& mat )
{
    auto [P, L, U, Q] = PLUQ( mat );
    auto check = P * L * U * Q;
    std::cout << "PLUQ check: \n";
    print( check );
    std::cout << '\n';
}



int main()
{
    using T = boost::multiprecision::mpfr_float_1000;
//    using T = double;

    using namespace boost::multiprecision;
//
//    Mat<float, 4, 4> mat;
//    mat( 0, 0 ) = 1;
//    mat( 0, 1 ) = 2;
//    mat( 0, 2 ) = 3;
//    mat( 0, 3 ) = 4;
//
//    mat( 1, 0 ) = 1;
//    mat( 1, 1 ) = 1;
//    mat( 1, 2 ) = 1;
//    mat( 1, 3 ) = 1;
//
//    mat( 2, 0 ) = 4;
//    mat( 2, 1 ) = 2;
//    mat( 2, 2 ) = 3;
//    mat( 2, 3 ) = 1;
//
//    mat( 3, 0 ) = 2;
//    mat( 3, 1 ) = 3;
//    mat( 3, 2 ) = 1;
//    mat( 3, 3 ) = 1;

    Mat<T, 3, 3> mat;
    mat( 0, 0 ) = 1;
    mat( 0, 1 ) = pow( T( 2 ), T( 20 ) );
    mat( 0, 2 ) = pow( T( 2 ), T( 40 ) );

    mat( 1, 0 ) = 2;
    mat( 1, 1 ) = pow( T( 2 ), T( 40 ) );
    mat( 1, 2 ) = pow( T( 2 ), T( 150 ) //    auto [P, L, U, Q] = PLUQ( mat );
//    print( P );
//    std::cout << '\n';
//    print( L );
//    std::cout << '\n';
//    print( U );
//    std::cout << '\n';
//    print( Q );
//    std::cout << '\n';
//
//    auto check = P * L * U * Q;
//    print( check );
//    std::cout << '\n';);

    mat( 2, 0 ) = pow( T( 2 ), T( 30 ) );
    mat( 2, 1 ) = pow( T( 2 ), T( 54 ) );
    mat( 2, 2 ) = pow( T( 2 ), T( 10 ) );

    print( mat );
    std::cout << '\n';

    std::cout << "Floating point:\n";
    tryLU( mat.cast<float>() );
    tryPLUQ( mat.cast<float>() );

    std::cout << "Multiprecision:\n";
    tryPLUQ( mat );

    return 0;
}
