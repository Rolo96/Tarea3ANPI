//
// Created by bryan on 08/03/18.
//

#include <functional>
#include <cmath>
namespace anpi {
    template <typename T>
    class encapsulator{
    private:
        static int f1count;
        static int f2count;
        static int f3count;
    public:
        encapsulator() {}
        int getF1count() const { return f1count; }
        int getF2count() const { return f2count; }
        int getF3count() const { return f3count; }
        static void setF1count(){f1count++;}
        static void setF2count(){encapsulator::f2count++;}
        static void setF3count(){encapsulator::f3count++;}

        static inline T sqr(const T x) { return x * x; }

        static T funct1(const T x) {
            setF1count();
            return std::abs(x) - std::exp(-x);
        }

        static T funct2(const T x) {
            setF2count();
            return std::exp(-x * x) - std::exp(-sqr(x - T(3)) / T(3));
        }

        static T funct3(const T x) {
            setF3count();
            return x * x - std::atan(x);
        }
    };

    template <typename T>
    int encapsulator<T>::f1count = 0;
    template <typename T>
    int encapsulator<T>::f2count = 0;
    template <typename T>
    int encapsulator<T>::f3count = 0;
}
