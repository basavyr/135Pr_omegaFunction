#include <iostream>
#include <fstream>
#include "../include/omega.h"

void newline()
{
    std::cout << "\n";
}

void app()
{
    std::cout << "App runs ok";
    newline();
}

template <typename T>
void printArrayToFile(std::ofstream &outstream, std::vector<T> &v)
{
    if (v.size())
    {
        outstream << "l1= {";
        for (auto id = 0; id < v.size(); ++id)
        {
            if (id == v.size() - 1)
            {
                outstream << v.at(id) << " };";
            }
            else
            {
                outstream << v.at(id) << " , ";
            }
        }
    }
    else
    {
        outstream << "l1= { 0};";
        outstream << std::endl;
    }
}

void generateOmegas(WobblingFrequency::X_Set &params, std::vector<double> &omegas)
{
    if (!omegas.size())
    {
        //only work if the entire array has real, valid numbers (physical solutions)
        int safety = 1;
        for (int theta = -180; theta <= 180 && safety; ++theta)
        {
            auto omega = WobblingFrequency::Omega(params, static_cast<double>(theta));
            if (isnan(omega) || omega == WobblingFrequency::ERROR)
            {
                safety = 0;
                return;
            }
            else
            {
                omegas.emplace_back(omega);
            }
        }
    }
    else
    {
        return;
    }
}

int main()
{
    WobblingFrequency::X_Set paramSet(20, 100, 40, 45.0 / 2.0, 11.0 / 2.0);
    WobblingFrequency::X_Set::printX_Set(paramSet);
    std::vector<double> omegas;
    generateOmegas(paramSet, omegas);

    std::ofstream gout;
    gout.open("../sources/omegas.dat", std::ios::trunc);

    printArrayToFile(gout, omegas);
    std::cout << WobblingFrequency::Omega(paramSet, 0.0);
}