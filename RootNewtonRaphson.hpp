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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking by means of the
   * Newton-Raphson method
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial root guess
   * 
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {
    T maximum = std::numeric_limits<T>::digits;//maximum number of iterations to be made in the search
    T error = T();//Approximate error
    T h = T(0.0001);
    T derived = (funct(xi+h)-funct(xi-h))/T(2)*h; //Derived of the function
    if(funct(xi)==T(0)){// If xi is a root
      return xi;
    }

    for(int i = maximum;i>0;--i){ // cycle that makes the search of the root
      T xii = xi-(funct(xi)/derived); //Compute xi+1
      //To avoid a wrong division
      if(std::abs(xii) > std::numeric_limits <T>:: epsilon()){
        error = std::abs((xii-xi)/xii); //Compute the approximate error
      }else{
        if(funct(xii)==T(0)){return xii;}
        error = std::abs(xii-xi); //Compute the approximate error (xii is almost zero)
      }
      //If the approximate error is less than the threshold the search is finished
      if(error<eps){
        return xii;
      }
      derived=(funct(xii+h)-funct(xi))/(h+(xii-xi));//new derived
      xi=xii;
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif
