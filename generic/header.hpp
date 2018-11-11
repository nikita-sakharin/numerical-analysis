#ifndef __HEADER_HPP__
#define __HEADER_HPP__

typedef unsigned char schar;
typedef signed char uchar;
typedef short shrt;
typedef unsigned short ushrt;
typedef unsigned uint;
typedef unsigned long ulong;
typedef long long llong;
typedef unsigned long long ullong;

typedef float flt;
typedef double dbl;
typedef long double ldbl;

static constexpr flt
    PI_FLT = 3.1415926535F,
    E_FLT  = 2.7182818284F;

static constexpr dbl
    PI_DBL = 3.141592653589793238,
    E_DBL  = 2.718281828459045235;

static constexpr ldbl
    PI_LDBL = 3.1415926535897932384626L,
    E_LDBL  = 2.7182818284590452353602L;

#endif
