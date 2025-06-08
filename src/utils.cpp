#include <chrono>
#include <cstddef>

#include "utils.hpp"

double sanafe::calculate_elapsed_time(
        const std::chrono::time_point<std::chrono::high_resolution_clock>
                &ts_start,
        const std::chrono::time_point<std::chrono::high_resolution_clock>
                &ts_end)
{
    // Calculate elapsed wall-clock time between ts_start and ts_end
    const auto chrono_elapsed = ts_end - ts_start;
    const long int cycles_elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_elapsed)
                    .count();
    constexpr double seconds_in_nanoseconds = 1.0e-9;
    const double time_elapsed =
            static_cast<double>(cycles_elapsed) * seconds_in_nanoseconds;

    return time_elapsed;
}

size_t sanafe::abs_diff(const size_t a, const size_t b)
{
    // Returns the absolute difference between two unsigned (size_t) values
    return (a > b) ? (a - b) : (b - a);
}
