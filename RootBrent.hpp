/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "iostream"
#include "Exception.hpp"

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the Brent's method.
   *
   * @param funct a std::function of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps) {
    if(xl>=xu){throw Exception("interval is reversed");}
    if(xl*xu>0){throw Exception("interval doesn't enclose the root");}
    T es =((xu-xl)/3)*eps;// new threshold
    T maximum = std::numeric_limits<T>::digits;//maximum number of iterations to be made in the search
    T fxl = funct(xl);
    T fxu = funct(xu);   
    T fxr = T();      // initialize
    if (std::abs(fxl) < std::abs(xu)){
      std::swap(xl,xu);
      std::swap(fxl,fxu);
    }
    T c = xl;// c now equals the largest magnitude of the lower and upper bounds
    T fc = fxl;// precompute function evalutation for point c by assigning it the same value as fxl
    bool mflag = true; // boolean flag used to evaluate in a if later
    T xr = T();// root variable
    T unsetflag = T();// Only used if mflag is unset (mflag == false)
    for (int i = 1; i < maximum; ++i){
      //If the approximate error is less than the threshold the search is finished
      if (std::abs(xu-xl) < es) {
        return xr;
      }
      if (fxl != fc && fxu != fc) { //the next approximation is compute with quadratic inverse interpolation
        xr=(xl*fxu*fc/((fxl-fxu)*(fxl-fc)))
              +(xu*fxl*fc/((fxu-fxl)*(fxu-fc)))
              +(c*fxl*fxu/((fc-fxl)*(fc-fxu)));
      }else{ //the next approximation is compute with secant
        xr=xu-fxu*(xu-xl)/(fxu-fxl);
      }
      //Condition to chose if use the bisection approximation
      if (( (xr<(3*xl+xu)*0.25)||(xr>xu))||(mflag&&(std::abs(xr-xu)>=(std::abs(xu-c)*0.5)))||
              (!mflag&&(std::abs(xr-xu)>=(std::abs(c-unsetflag)*0.5)))||(mflag&&(std::abs(xu-c)<es))||
              (!mflag&&(std::abs(c-unsetflag)<es))){
        xr = (xl+xu)*0.5;
        mflag = true;
      }else{
        mflag = false;
      }
      fxr = funct(xr);
      unsetflag = c;// first time unsetflag is being used (wasnt used on first iteration because mflag was set)
      c = xu;// set c equal to upper bound
      fc = fxu;
      if(fxl*fxr < 0){// if fa and fs have opposite signs so:
        xu = xr;
        fxu = fxr;
      }else{
        xl = xr;
        fxl = fxr;
      }
      if (std::abs(fxl) < std::abs(fxu)){ // if magnitude of fa is less than magnitude of fb
        std::swap(xl,xu);
        std::swap(fxl,fxu);
      }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }
}
  
#endif

