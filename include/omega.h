#ifndef OMEGA_HH
#define OMEGA_HH

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

class WobblingFrequency
{
public:
    static constexpr double PI = 3.14159265358979;
    static constexpr int ERROR = 6969;

public:
    //define the set of parameters which are predefined for calculating omega and omega prime
    struct X_Set
    {
        double I;
        double j;
        double I1;
        double I2;
        double I3;
        X_Set(double i1, double i2, double i3, double spin, double oddspin)
        {
            I = static_cast<double>(spin);
            j = static_cast<double>(oddspin);
            I1 = static_cast<double>(i1);
            I2 = static_cast<double>(i2);
            I3 = static_cast<double>(i3);
        }
        static void printX_Set(X_Set &params)
        {
            auto A1 = inertiaFactor(params.I1);
            auto A2 = inertiaFactor(params.I2);
            auto A3 = inertiaFactor(params.I3);
            std::cout << "I= " << params.I << "\n";
            std::cout << "j= " << params.j << "\n";
            std::cout << "A1= " << A1 << "\n";
            std::cout << "A2= " << A2 << "\n";
            std::cout << "A3= " << A3 << "\n";
        }
    };
    //define helper functions which are using inside the omega defintion and parameters
    static double inertiaFactor(double Ik);
    static double j_Component(int k, double theta, const double j);
    //Omega is a function which depends on several parameters (see @README for information about the parameters)
    //Omega represents the wobbling frequency associated with the energy spectrum of the Rotor given by the total Hamiltonian H_rot
    static double Omega(X_Set &paramSet, double theta);
    //Omega' is the oscillator frequency associated with a triaxial rotor which is moving in a potential, only around the minimum point.
    static double OmegaPrime(X_Set &paramSet, double theta);
};

#endif // OMEGA_HH
