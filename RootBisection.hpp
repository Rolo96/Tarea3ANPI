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
#include <iostream>
#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the bisection method.
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
  T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

      if (xl>=xu){throw Exception("inteval is reversed");}

      T maximum = std::numeric_limits<T>::digits;//maximum number of iterations to be made in the search
      T xr = xl; //Start with xr = xl
      T fl = funct(xl); //The function is evaluated in xl
      T error = T(100);//Approximate error

      T es = ((xu-xl)/2)*eps;

      // cycle that does the searches dividing the interval in half in each iteration
      for(int i = maximum;i>0;--i){
          T lastXR(xr);
          xr = (xu+xl)/T(2); // new xr
          T fr = funct(xr); //he function is evaluated in the new xl
          error = std::abs(xr-lastXR);// The approximate error is calculated

          T decision = fl*fr; //  know in which half of the interval to keep looking for the root

          if(decision<T(0)){
              xu=xr; // There is change of sign in xl
          }else if(decision>T(0)){
              xl=xr; // There is no change of sign in xl
              fl=fr;
          }else{ // fl o fr is zero, that is, xl or xr is a root
              error=T(0);

              xr = (std::abs(fl) < eps)?xl:xr; // if fl is zero xr=xl, else xr=xr
          }
          if(error<es){ // If the approximate error is less than es
              if (funct(xr)<=eps) return xr; // xr is a root
              else throw Exception("unenclosed root");
          }
      }
      // Return NaN if no root was found
      throw Exception("unenclosed root");
      //return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif