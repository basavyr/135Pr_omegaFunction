#include <iostream>
#include <fstream>
#include "../include/omega.h"
#include <string>

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
void printArrayToFile(std::ofstream &outstream, std::vector<T> &v, const std::string &plotName)
{
    if (v.size())
    {
        outstream << plotName;
        for (auto id = 0; id < v.size(); ++id)
        {
            if (id == v.size() - 1)
            {
                outstream << v.at(id) << " };";
                outstream << "\n";
            }
            else
            {
                outstream << v.at(id) << " , ";
            }
        }
    }
    else
    {
        std::string localplot = plotName;
        localplot.append("0 };");
        outstream << localplot;
        outstream << "\n";
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
                std::cout << omega << "\n";
                std::cout << "Parameters are not providing physical solutions to the wobbling frequencies"
                          << "\n";
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

void generateOmegaPrimes(WobblingFrequency::X_Set &params, std::vector<double> &omegas)
{
    if (!omegas.size())
    {
        //only work if the entire array has real, valid numbers (physical solutions)
        int safety = 1;
        for (int theta = -180; theta <= 180 && safety; ++theta)
        {
            auto omega = WobblingFrequency::OmegaPrime(params, static_cast<double>(theta));
            if (isnan(omega) || omega == WobblingFrequency::ERROR)
            {
                std::cout << omega << "\n";
                std::cout << "Parameters are not providing physical solutions to the wobbling frequencies"
                          << "\n";
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

void generateOmegaContainers(std::vector<double> &spins)
{
    auto I1 = 20;
    auto I2 = 100;
    auto I3 = 40;
    auto j = 5.5;

    std::string ListPlot = "ListPlot[{";

    std::ofstream gout;
    gout.open("../sources/omegas.dat", std::ios::trunc);

    int plotIdx = 1;
    for (auto id = 0; id < spins.size(); ++id)
    {
        auto I = spins.at(id);
        WobblingFrequency::X_Set paramSet(I1, I2, I3, I, j);
        std::vector<double> omegas;
        generateOmegas(paramSet, omegas);
        gout << "l" << plotIdx << "= { ";
        for (auto omegaIdx = 0; omegaIdx < omegas.size(); ++omegaIdx)
        {
            if (omegaIdx == omegas.size() - 1)
            {
                gout << omegas.at(omegaIdx) << " };";
                gout << "\n";
            }
            else
            {
                gout << omegas.at(omegaIdx) << " , ";
            }
        }
        if (plotIdx == spins.size())
        {
            //construct the name of the current plot based on the idx
            std::string currentPlot = "l";
            currentPlot.append(std::to_string(plotIdx));
            ListPlot.append(currentPlot);
            ListPlot.append("},Joined->True]");
        }
        else
        {
            std::string currentPlot = "l";
            currentPlot.append(std::to_string(plotIdx));
            ListPlot.append(currentPlot);
            ListPlot.append(",");
        }
        plotIdx++;
    }
    gout << ListPlot;
    gout << "\n";
}

void generateOmegaComparison(double spin)
{
    auto I1 = 20;
    auto I2 = 100;
    auto I3 = 40;
    auto j = 6.5;
    WobblingFrequency::X_Set paramSet(I1, I2, I3, spin, j);
    std::vector<double> omegas;
    std::vector<double> omegaprimes;
    for (int theta = 0; theta <= 180; ++theta)
    {
        auto currentOmega = WobblingFrequency::Omega(paramSet, static_cast<double>(theta));
        auto currentOmegaPrime = WobblingFrequency::OmegaPrime(paramSet, static_cast<double>(theta));
        if (isnan(currentOmega) || isnan(currentOmegaPrime) || currentOmega == WobblingFrequency::ERROR || currentOmegaPrime == WobblingFrequency::ERROR)
        {
            std::cout << "Comparison is invalid due to non-physical solutions in the wobbling frequencies"
                      << "\n";
            theta = 200;
            return;
        }
        omegas.emplace_back(currentOmega);
        omegaprimes.emplace_back(currentOmegaPrime);

        std::string omegaPlot = "omega={ ";
        std::string omegaPrimePlot = "omegaPrime={ ";

        std::ofstream gout;
        gout.open("../sources/omegas.dat", std::ios::trunc);

        printArrayToFile(gout, omegas, omegaPlot);
        printArrayToFile(gout, omegaprimes, omegaPrimePlot);
    }
}

void generateOmegasFromSpin(std::ofstream &out, double spin)
{
    //generate the omega containers from a predefined set of parameters

    //the parameters obtained from the C++ best fit results
    WobblingFrequency::X_Set X_CPP(89.0, 12.0, 48.0, spin, 5.5);

    //the parameters obtained from the Mathematica NSolve system of equations
    WobblingFrequency::X_Set X_Math(13.52, 101.759, 52.936, spin, 5.5);

    std::vector<double> omega_CPP;
    std::vector<double> omegaPrime_CPP;

    std::vector<double> omega_Math;
    std::vector<double> omegaPrime_Math;

    generateOmegas(X_CPP, omega_CPP);
    generateOmegaPrimes(X_CPP, omegaPrime_CPP);

    generateOmegas(X_Math, omega_Math);
    generateOmegaPrimes(X_Math, omegaPrime_Math);

    std::string cppPlot = "omegaCPP= {";
    std::string cppPlotPrime = "omegaCPPPrime= {";

    std::string mathPlot = "omegaMath= {";
    std::string mathPlotPrime = "omegaMathPrime= {";

    // printArrayToFile(out, omega_Math, mathPlot);
    // printArrayToFile(out, omegaPrime_Math, mathPlotPrime);

    printArrayToFile(out, omega_CPP,cppPlot);
    printArrayToFile(out, omegaPrime_CPP, cppPlotPrime);
}

int main()
{
    std::vector<double> spins = {11.5, 15.5, 19.5, 23.5};
    // generateOmegaContainers(spins);
    // generateOmegaComparison(22.5);

    std::ofstream out;
    out.open("../sources/omegas.dat", std::ios::trunc);
    // generateOmegasFromSpin(out, 22.5);
    generateOmegasFromSpin(out, 19.0 / 2.0);
}