#pragma once
/** 
 * spline_interpolator.h - Utility classes for 1D cubic spline interpolation
 *
 * Copyright (C) 2017 University of Oxford 
 */

/*  CCOPYRIGHT */

#include <vector>
#include <string>

/** 
 * Parameters for cubic spline starting at x=x0
 * 
 * Equation is d(x-x0)^3 + c(x-x0)^2 + b(x-x0) + 1
 */
struct Spline {
    Spline(double a=0, double b=0, double c=0, double d=0, double x0=0)
      : a(a), b(b), c(c), d(d), x0(x0) {}

    /** Evaluate the equation */
    double operator()(double x) const;

    /** Return human-readable equation to plug into FooPlot etc. */
    std::string eqn() const;

    double a, b, c, d, x0;
};

/**
 * Base class for anything which does spline interpolation
 */
class SplineInterpolator
{
public:
    double operator()(double x) const;
protected:
    std::vector<Spline> m_splines;
};

/**
 * PCHIP spline interpolator, as used in Matlab version of the model
 */
class PchipInterpolator : public SplineInterpolator
{
public:
    PchipInterpolator(std::vector<double> &x, std::vector<double> &y);
private:
    Spline get_spline(double x1, double y1, double m1, double x2, double y2, double m2);
    double edge_case(double h0, double h1, double m0, double m1);
    std::vector<double> get_derivatives(std::vector<double> &x, std::vector<double> &y);
};

/**
 * 'Natural' spline interpolator for comparison
 */
class NaturalSplineInterpolator : public SplineInterpolator
{
public:
    NaturalSplineInterpolator(std::vector<double> &x, std::vector<double> &y);
};
