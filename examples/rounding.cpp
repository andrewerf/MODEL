//
// Created by Andrey Aralov on 9/23/24.
//
#include <iostream>
#include <cmath>
#include <iomanip>

int main()
{
    double a = 1;
    double b = pow( double( 2 ), 90 );
    double c = a + b;

    std::cout << std::setprecision( 10000 ) << a << '\n' << b << '\n' << c << '\n';
}