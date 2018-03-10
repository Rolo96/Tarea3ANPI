/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>
#include "benchmarkRootFinders.cpp"

#include "RootInterpolation.hpp"
#include "RootBrent.hpp"
#include "RootNewtonRaphson.hpp"
#include "RootBisection.hpp"
#include "RootSecant.hpp"

#include "vector"
#include "PlotPy.hpp"

template <typename T>
void solve(const std::function<T(const std::function<T(T)>&,
                                T,
                                T,
                                const T)>& solver, std::vector<T> &counts1, std::vector<T> &counts2, std::vector<T> &counts3){
    anpi::encapsulator<T> e;
    int index = 0;
    for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {
        e.setF1count(-e.getF1count());
        solver(e.template funct1,T(0),T(2),eps);
        counts1[index] = T(e.getF1count());

        e.setF2count(-e.getF2count());
        solver(e.template funct2,T(0),T(2),eps);
        counts2[index] = T(e.getF2count());

        e.setF3count(-e.getF3count());
        solver(e.template funct1,T(0),T(2),eps);
        counts3[index++] = T(e.getF3count());
    }

}

template <typename T>
void solve2(const std::function<T(const std::function<T(T)>&,
                                 T,
                                 const T)>& solver, std::vector<T> &counts1, std::vector<T> &counts2, std::vector<T> &counts3){
    anpi::encapsulator<T> e;
    int index  = 0;
    for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(10)) {
        e.setF1count(-e.getF1count());
        solver(e.template funct1,T(0),eps);
        counts1[index] = T(e.getF1count());

        e.setF2count(-e.getF2count());
        solver(e.template funct2,T(2),eps);
        counts2[index] = T(e.getF2count());

        e.setF3count(-e.getF3count());
        solver(e.template funct3,T(0),eps);
        counts3[index++] = T(e.getF3count());
    }

}

int main() {
    anpi::Plot2d<float> *p1 = new anpi::Plot2d<float>();
    p1->initialize(1);

    anpi::Plot2d<double> *p2d = new anpi::Plot2d<double>();
    p2d->initialize(1);

    float epsf = 0.1f;
    double epsd = 0.1;

    std::vector<float> Epsf(6);
    std::vector<double > Epsd(7);
    for (int i = 0; i < 6; ++i) { //esto porque los vectores del epsilon son iguales
                                    // para todos entonces se hace solo uno para float y double
                                        // No se si sera buena idea
        Epsd[i] = epsd;
        Epsf[i] = epsf;
        epsd /= 10.0;
        epsf /= 10.0f;
    }
    epsd /= 10.0;
    Epsd[6] = epsd;

    std::vector<float> countsf1(6);
    std::vector<float> countsf2(6);
    std::vector<float> countsf3(6);
    std::vector<double> countsd1(7);
    std::vector<double> countsd2(7);
    std::vector<double> countsd3(7);

    solve<float>(anpi::rootInterpolation<float>,countsf1,countsf2,countsf3);
    p1->figure(1);
    p1->plot(countsf1,Epsf,"Interpolacion","blue");
    p1->figure(2);
    p1->plot(countsf2,Epsf,"Interpolacion","blue");
    p1->figure(3);
    p1->plot(countsf3,Epsf,"Interpolacion","blue");
    solve<float>(anpi::rootBisection<float>,countsf1,countsf2,countsf3);
    p1->figure(1);
    p1->plot(countsf1,Epsf,"Biseccion","red");
    p1->figure(2);
    p1->plot(countsf2,Epsf,"Biseccion","red");
    p1->figure(3);
    p1->plot(countsf3,Epsf,"Biseccion","red");
    solve<float>(anpi::rootSecant<float>,countsf1,countsf2,countsf3);
    p1->figure(1);
    p1->plot(countsf1,Epsf,"Secante","yellow");
    p1->figure(2);
    p1->plot(countsf2,Epsf,"Secante","yellow");
    p1->figure(3);
    p1->plot(countsf3,Epsf,"Secante","yellow");
    solve<float>(anpi::rootBrent<float>,countsf1,countsf2,countsf3);
    p1->figure(1);
    p1->plot(countsf1,Epsf,"Brent","green");
    p1->figure(2);
    p1->plot(countsf2,Epsf,"Brent","green");
    p1->figure(3);
    p1->plot(countsf3,Epsf,"Brent","green");
    solve2<float>(anpi::rootNewtonRaphson<float>,countsf1,countsf2,countsf3);
    p1->figure(1);
    p1->plot(countsf1,Epsf,"Newton","orange");
    p1->figure(2);
    p1->plot(countsf2,Epsf,"Newton","orange");
    p1->figure(3);
    p1->plot(countsf3,Epsf,"Newton","orange");

    solve<double>(anpi::rootInterpolation<double>,countsd1,countsd2,countsd3);
    p2d->figure(4);
    p2d->plot(countsd1,Epsd,"Interpolacion","blue");
    p2d->figure(5);
    p2d->plot(countsd2,Epsd,"Interpolacion","blue");
    p2d->figure(6);
    p2d->plot(countsd3,Epsd,"Interpolacion","blue");
    solve<double>(anpi::rootBisection<double>,countsd1,countsd2,countsd3);
    p2d->figure(4);
    p2d->plot(countsd1,Epsd,"Biseccion","red");
    p2d->figure(5);
    p2d->plot(countsd2,Epsd,"Biseccion","red");
    p2d->figure(6);
    p2d->plot(countsd3,Epsd,"Biseccion","red");
    solve<double>(anpi::rootSecant<double>,countsd1,countsd2,countsd3);
    p2d->figure(4);
    p2d->plot(countsd1,Epsd,"Secante","yellow");
    p2d->figure(5);
    p2d->plot(countsd2,Epsd,"Secante","yellow");
    p2d->figure(6);
    p2d->plot(countsd3,Epsd,"Secante","yellow");
    solve<double>(anpi::rootBrent<double>,countsd1,countsd2,countsd3);
    p2d->figure(4);
    p2d->plot(countsd1,Epsd,"Brent","green");
    p2d->figure(5);
    p2d->plot(countsd2,Epsd,"Brent","green");
    p2d->figure(6);
    p2d->plot(countsd3,Epsd,"Brent","green");
    solve2<double>(anpi::rootNewtonRaphson<double>,countsd1,countsd2,countsd3);
    p2d->figure(4);
    p2d->plot(countsd1,Epsd,"Newton","orange");
    p2d->figure(5);
    p2d->plot(countsd2,Epsd,"Newton","orange");
    p2d->figure(6);
    p2d->plot(countsd3,Epsd,"Newton","orange");

    p1->show();

    return EXIT_FAILURE;
}
  
  
