#include <random>
#include <util.h>

namespace Random {
    extern std::mt19937 mt { std::random_device {}() };
    double get(double min, double max) {
        std::uniform_real_distribution dis { min, max };
        return dis(mt);
    }
}