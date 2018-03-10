//
// Created by bryan on 08/03/18.
//

#include <functional>
#include <cmath>
namespace anpi {
    template <typename T>
    /**
     * Encapsulator class
     */
    class encapsulator{
    private:
        //Function Counters
        static int f1count;
        static int f2count;
        static int f3count;
    public:
        // Class Constructor
        encapsulator() {}
        //Counter Getters and Setters
        int getF1count() const { return f1count; }
        int getF2count() const { return f2count; }
        int getF3count() const { return f3count; }
        static void setF1count(int num){f1count+=num;}
        static void setF2count(int num){f2count+=num;}
        static void setF3count(int num){f3count+=num;}
        //Encapsulated functions
        static inline T sqr(const T x) { return x * x; }
        static T funct1(const T x){setF1count(1);return std::abs(x) - std::exp(-x);}
        static T funct2(const T x){setF2count(1);return std::exp(-x * x) - std::exp(-sqr(x - T(3)) / T(3));}
        static T funct3(const T x){setF3count(1);return x * x - std::atan(x);}
    };
    //Counters Initialization
    template <typename T>
    int encapsulator<T>::f1count = 0;
    template <typename T>
    int encapsulator<T>::f2count = 0;
    template <typename T>
    int encapsulator<T>::f3count = 0;
}
