#ifndef SLOC_UTILS_H
#define SLOC_UTILS_H

#define _PRINT_VALUE(x) \
    do { std::cout << #x << " = " << (x) << std::endl; } while (false)

#define _PRINT_EXEC(x) \
    do { std::cout << #x << std::endl; x; } while (false)

bool is_nan(double x);

#endif
