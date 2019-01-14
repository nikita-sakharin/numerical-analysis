#ifndef __HEADER_HPP__
#define __HEADER_HPP__

typedef unsigned char uchar;
typedef signed char schar;
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
    E_FLT    = 2.7182818284F,
    PI_FLT   = 3.1415926535F,
    PI_2_FLT = 1.5707963267F,
    LN2_FLT  = 6.9314718055E-1F;

static constexpr dbl
    E_DBL    = 2.7182818284590452353,
    PI_DBL   = 3.1415926535897932384,
    PI_2_DBL = 1.5707963267948966192,
    LN2_DBL  = 6.9314718055994530941E-1;

static constexpr ldbl
    E_LDBL    = 2.7182818284590452353602L,
    PI_LDBL   = 3.1415926535897932384626L,
    PI_2_LDBL = 1.5707963267948966192313L,
    LN2_LDBL  = 6.9314718055994530941723E-1L;

#endif
