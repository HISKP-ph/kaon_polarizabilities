#ifndef TYPE_ALIASES_H
#define TYPE_ALIASES_H

#include <complex>
#include <functional>

///Define type aliases for easier access.
namespace type_aliases {
/// define Complex double
using Complex = std::complex<double>;
/// Complex function with real arguments
using Curve = std::function<Complex(double)>;
/// Complex function with complex aruments
using CFunction = std::function<Complex(Complex)>;
} // type_aliases

#endif // TYPE_ALIASES_H
