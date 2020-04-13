#include "../include/omega.h"

double WobblingFrequency::inertiaFactor(double Ik)
{
    if (Ik)
        return static_cast<double>(1.0 / (2.0 * Ik));
    return 0;
}

double WobblingFrequency::j_Component(int k, double theta, const double j)
{
    auto theta_Degrees = theta * PI / 180.0;
    if (k == 1)
        return j * cos(theta_Degrees);
    return j * sin(theta_Degrees);
}

double WobblingFrequency::Omega(WobblingFrequency::X_Set &paramSet, double theta)
{
    auto A1 = inertiaFactor(paramSet.I1);
    auto A2 = inertiaFactor(paramSet.I2);
    auto A3 = inertiaFactor(paramSet.I3);

    //stop if the moments of inertia are null
    if (!A1 && !A2 && !A3)
        return ERROR;

    auto j1 = j_Component(1, theta, paramSet.j);
    auto j2 = j_Component(2, theta, paramSet.j);
    auto I = paramSet.I;

    //construct the terms of the frequency function
    auto term1 = (2.0 * I + 1.0) * (A2 - A1 - (A2 * j2) / I) - 2.0 * A1 * j1;
    if (isnan(term1))
        return ERROR;
    auto term2 = (2.0 * I + 1.0) * (A3 - A1) - 2.0 * A1 * j1;
    if (isnan(term2))
        return ERROR;
    auto term3 = (A3 - A1) * (A2 - A1 - (A2 * j2) / I);
    if (isnan(term3))
        return ERROR;

    //construct the final form of the wobbling frequency
    auto omega = static_cast<double>(sqrt(term1 * term2 - term3));
    if (!isnan(omega))
        return omega;
    return ERROR;
}

double WobblingFrequency::OmegaPrime(WobblingFrequency::X_Set &paramSet, double theta)
{
    auto A1 = inertiaFactor(paramSet.I1);
    auto A2 = inertiaFactor(paramSet.I2);
    auto A3 = inertiaFactor(paramSet.I3);

    //stop if the moments of inertia are null
    if (!A1 && !A2 && !A3)
        return ERROR;

    auto j1 = j_Component(1, theta, paramSet.j);
    auto j2 = j_Component(2, theta, paramSet.j);
    auto I = paramSet.I;

    //construct the terms of the frequency function
    auto term1 = (2.0 * I + 1.0) * (A2 - A1 - (A2 * j2) / I) + 2 * A1 * j1;
    auto term2 = (2.0 * I + 1.0) * (A3 - A1) + 2.0 * A1 * j1;
    auto term3 = (A3 - A1) * (A2 - A1 - (A2 * j2) / I);

    //construct the final form of the wobbling frequency
    auto omega = static_cast<double>(sqrt(term1 * term2 - term3));
    if (!isnan(omega))
        return omega;
    return ERROR;
}