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

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {
    const int maximum = std::numeric_limits<T>::digits;
    T error = T();
    T root_i = xi;
    T root_ii = xii;

    for(int i=maximum; i>0; --i){

      T y_i = funct(root_i);//evaluation of the function in Xi-1
      T y_ii = funct(root_ii);//evaluation of the function in Xi

      T root_iii = root_ii - (y_ii*(root_i-root_ii))/ (y_i-y_ii);//next iteration

      if(std::abs(root_iii) > std::numeric_limits<T>::epsilon()){//avoid division by zero
        error = std::abs((root_iii-root_ii)/root_iii) * T(100);
      }

      if(error < eps){return root_iii;} //stop condition
      root_i = root_ii;// Xi-1 = Xi
      root_ii = root_iii;//Xi=Xi+1
    }
    return std::numeric_limits<T>::quiet_NaN();
  }

}
#endif