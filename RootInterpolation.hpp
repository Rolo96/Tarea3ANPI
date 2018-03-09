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

#include "Exception.hpp"

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], by means of the interpolation method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootInterpolation(const std::function<T(T)>& funct,T xl,T xu,const T eps) {
    if(xl>=xu){throw Exception("interval is reversed");}
    if(funt(xl)*funt(xu)>0){throw Exception("interval doesn't enclose the root");}
    using std::abs;
    T maximum = std::numeric_limits<T>::digits;//maximum number of iterations to be made in the search
    T error=T();//Approximate error
    T xr = xl; //Start with xr = xl
    //the extremes are evaluated
    T fl = funct(xl);
    T fu = funct(xu);

    //To check that the same extreme is not chosen more than twice
    int upStuck = 0;
    int downStuck = 0;

    for (int i=maximum; i>0;--i){ //cycle that makes the search of the root
      T lastXr = xr; //Save the previous approximation
      xr = xu-(fu*(xl-xu))/(fl-fu);//compute the actual approximation
      T fr=funct(xr);//Evaluate the actual approximation

      //To avoid a wrong division
      if (abs(xr)>std::numeric_limits<T>::epsilon()){
        error = std::abs((xr-lastXr)/xr)*T(100);
      }

      T condition = fl*fr;//Condition to chose the extreme

      //Is nearer to lower extreme
      if (condition<T(0)){
        xu=xr;
        fu=fr;
        upStuck=0;
        downStuck++;
        if(downStuck>=2){fl/=T(2);}
      }
      //Is nearer to upper extreme
      else if (condition>T(0)){
        xl=xr;
        fl=fr;
        downStuck=0;
        upStuck++;
        if(upStuck>=2){fu/=2;}
      }
        //Bingo! A root
      else{
        error=T(0);//Error is zero
      }
      //If the approximate error is less than the threshold the search is finished
      if (error<eps)return xr;
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

