//
// Created by Andrey Aralov on 10/7/24.
//
#include "QR.hpp"
#include "io.hpp"


int main()
{
    Mat<double, 4, 2> mat;
    mat( 0, 0 ) = -7;
    mat( 0, 1 ) = 21;
    mat( 1, 0 ) = -4;
    mat( 1, 1 ) = 26;
    mat( 2, 0 ) = -4;
    mat( 2, 1 ) = -2;
    mat( 3, 0 ) = 0;
    mat( 3, 1 ) = 7;

//    Mat<float, 3, 2> mat;
//    mat( 0, 0 ) = 3;
//    mat( 0, 1 ) = -3;
//    mat( 1, 0 ) = 4;
//    mat( 1, 1 ) = -4;
//    mat( 2, 0 ) = 0;
//    mat( 2, 1 ) = 40;

    auto [Q, R] = QR_Given( mat );
    print( Q );
    std::cout << '\n';
    print( R );
    std::cout << '\n';


    auto check = Q * R;
    print(check);
}