//
// Created by bryan on 08/03/18.
//

#include <functional>
#include <cmath>
namespace anpi {
    template <typename T>
    class encapsulator : std::function<T(T)>  {
    private:
        int f1count;
        int f2count;
        int f3count;
    public:
        encapsulator() {
            f1count=0;
            f2count=0;
            f3count=0;
        }
        int getF1count() const { return f1count; }
        int getF2count() const { return f2count; }
        int getF3count() const { return f3count; }

        inline T sqr(const T x) { return x * x; }

        T funct1(const T x) {
            f1count++;
            return std::abs(x) - std::exp(-x);
        }

        T funct2(const T x) {
            f2count++;
            return std::exp(-x * x) - std::exp(-sqr(x - T(3)) / T(3));
        }

        T funct3(const T x) {
            f3count++;
            return x * x - std::atan(x);
        }
    };
}
