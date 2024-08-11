/**  const1.h
 **  Copyright 2022-2024 Dmitry Boldyrev
 **  All rights reserved.
 **  
 **  This file is part of OFFT C++ class located at:
 **  https://github.com/mewza/offt/
 ** 
 **  GITHUB: https://github.com/mewza
 **  Email: subband@protonmail.com
 **
 **  A very useful set of tools for using SIMD on OS X, 
 **  but can probably easily port this to any clang++ or gcc compiler
 **
 **  LICENSE: Free for Commercial and Shareware use, however credit
 **  must be given to me in your app's about box.
 **/

#pragma once

#include <signal.h>
#include <malloc/malloc.h>
#include <simd/simd.h>
#include <type_traits>
#include <algorithm>
#include <map>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <functional>

#define MSSFLOAT_IS_FLOAT (sizeof(mssFloat) == sizeof(float))

template<typename T>
concept IsVector = requires { T{}[0]; };

#define NOT_VECTOR(Z) (std::is_same_v<Z, float> || std::is_same_v<Z, double>)
#define IS_VECTOR(Z)  (IsVector<Z>)


static inline constexpr int alignFor(int typeSize, int nelem) {
   // if (typeSize > 2)
     //   return 16;
    return typeSize * nelem;
}

template<typename ZZ>
using SimdBase = decltype([] {
    using std::is_same_v;

    if constexpr(
        is_same_v<ZZ, int> ||
        is_same_v<ZZ, float> ||
        is_same_v<ZZ, double>
    )
        return ZZ{};
    else
        return ZZ{}[0] + 0;
}());

template<typename ZZ> struct SimdInfo {
    using Base = SimdBase<ZZ>;
    static constexpr int size = sizeof(ZZ) / sizeof(Base);
};

template <class ZZ>
inline constexpr int SimdSize = SimdInfo<ZZ>::size;

template<typename Z, int size>
using Simd = Z __attribute__((ext_vector_type(size),aligned(alignFor(sizeof(Z),size))));

template<typename ZZ, typename NewBase>
using SimdSame = Simd<NewBase, SimdInfo<ZZ>::size>;

template<typename ZZ, typename NewBase>
using SimdSameHalf = Simd<NewBase, SimdInfo<ZZ>::size/2>;


#ifndef D_ZFLOAT
#define D_ZFLOAT

#define DECL_ZFLOAT(TYPE) \
typedef TYPE zfloat; \
typedef simd_##TYPE##8 zfloat8; \
typedef simd_##TYPE##4 zfloat4; \
typedef simd_##TYPE##2 zfloat2; \
static inline zfloat8 make_zfloat8(zfloat4& a, zfloat4& b) { return simd_make_##TYPE##8(a,b); } \
static inline zfloat4 make_zfloat4(zfloat2& a, zfloat2& b) { return simd_make_##TYPE##4(a,b); } \
static inline zfloat2 make_zfloat2(zfloat& a, zfloat& b) { return simd_make_##TYPE##2(a,b); }

DECL_ZFLOAT(double)

#endif


#define ZFLOAT_IS_FLOAT constexpr(sizeof(zfloat) == sizeof(float))
#define FUNDEMENTAL(T) std::is_same_v<T, float> || std::is_same_v<T, double>

typedef float mssFloat;

typedef simd_float8 mssFloat8;
typedef simd_float4 mssFloat4;
typedef simd_float2 mssFloat2;


typedef simd_float2 float2v;
typedef simd_float4 float4v;
typedef simd_float8 float8v;

typedef simd_double2 double2v;
typedef simd_double4 double4v;
typedef simd_double8 double8v;

typedef simd_long8 long8v;
typedef simd_long4 long4v;
typedef simd_long2 long2v;

typedef simd_int2 int2v;
typedef simd_int4 int4v;
typedef simd_int8 int8v;


template<class T, size_t N> struct ZArray;
template<class T, size_t N> struct ZArray {
    
    using IT = SimdSame<T,int>;
    using T8 = Simd<SimdBase<T>,8>;
    using T4 = Simd<SimdBase<T>,4>;
    using T2 = Simd<SimdBase<T>,2>;

public:
    T dd[N];
    
    class ZProxy
    {
    public:
        ZArray &a;
        IT ii;
        
        ZProxy(ZArray &a, IT& i) : a(a), ii(i) {}
        T& operator = (T& v) {
            if constexpr( SimdSize<T> == 8 ) {
                static T8 tmp;
                a.dd[ii[0]][0]=v[0]; a.dd[ii[1]][1]=v[1]; a.dd[ii[2]][2]=v[2]; a.dd[ii[3]][3]=v[3];
                a.dd[ii[4]][4]=v[4]; a.dd[ii[5]][5]=v[5]; a.dd[ii[6]][6]=v[6]; a.dd[ii[7]][7]=v[7];
                tmp = T{ a.dd[ii[0]][0], a.dd[ii[1]][1], a.dd[ii[2]][2], a.dd[ii[3]][3], a.dd[ii[4]][4], a.dd[ii[5]][5], a.dd[ii[6]][6], a.dd[ii[7]][7] };
                return tmp;
            } else if constexpr( SimdSize<T> == 4 ) {
                static T4 tmp;
                a.dd[ii[0]][0]=v[0]; a.dd[ii[1]][1]=v[1]; a.dd[ii[2]][2]=v[2]; a.dd[ii[3]][3]=v[3];
                tmp = T{ a.dd[ii[0]][0], a.dd[ii[1]][1], a.dd[ii[2]][2], a.dd[ii[3]][3] };
                return tmp;
            } else if constexpr( SimdSize<T> == 2 ) {
                static T2 tmp;
                a.dd[ii[0]][0]=v[0]; a.dd[ii[1]][1]=v[1];
                tmp = T{ a.dd[ii[0]][0], a.dd[ii[1]][1] };
                return tmp;
            }
        }
        
        inline operator T() const 
        {
            if constexpr( SimdSize<T> == 8 ) {
                return T{ a.dd[ii[0]][0], a.dd[ii[1]][1], a.dd[ii[2]][2], a.dd[ii[3]][3], a.dd[ii[4]][4], a.dd[ii[5]][5], a.dd[ii[6]][6], a.dd[ii[7]][7] };
            } else if constexpr( SimdSize<T> == 4 ) {
                return T{ a.dd[ii[0]][0], a.dd[ii[1]][1], a.dd[ii[2]][2], a.dd[ii[3]][3] };
            } else if constexpr( SimdSize<T> == 2 ) {
                return T{ a.dd[ii[0]][0], a.dd[ii[1]][1] };
            }
        }
    };
    ZProxy operator[] (IT& i) { return ZProxy(*this, i); }
    ZArray() {}
};

typedef std::conditional<sizeof(mssFloat) == sizeof(float), float2v, double2v>::type mssFloat2;
typedef std::conditional<sizeof(mssFloat) == sizeof(float), float4v, double4v>::type mssFloat4;
typedef std::conditional<sizeof(mssFloat) == sizeof(float), float8v, double8v>::type mssFloat8;

static __inline float8v make_mssFloat8(const float4v &a, const float4v &b) {
    return simd_make_float8(a,b);
}

static __inline double8v make_mssFloat8(const double4v &a, const double4v &b) {
    return simd_make_double8(a,b);
}

static __inline float4v make_mssFloat4(const float2v &a, const float2v &b) {
    return simd_make_float4(a,b);
}

static __inline double4v make_mssFloat4(const double2v &a, const double2v &b) {
    return simd_make_double4(a,b);
}

static __inline float2v make_mssFloat2(const float &a, const float &b) {
    return simd_make_float2(a,b);
}

static __inline double2v make_mssFloat2(const double &a, const double &b) {
    return simd_make_double2(a,b);
}

//typedef mssFloat cmplx[2];

#ifndef CMPLX_T_TYPE
#define CMPLX_T_TYPE

template <typename T> struct cmplxT;

#endif // CMPLX_T_TYPE

class const1 {
public:
    float v;
    constexpr const1 (mssFloat value) : v(value) {}
    inline operator mssFloat2() const {
        return mssFloat2{v, v};
    }
    inline operator mssFloat4() const {
        return mssFloat4{v, v, v, v};
    }
    inline operator mssFloat8() const {
        return mssFloat8{v, v, v, v, v, v, v, v};
    }
    inline operator float() const {
        return v;
    }
    inline constexpr const1 operator-() const { return const1(-v); }
};

#define shufflevector __builtin_shufflevector
#define convertvector __builtin_convertvector

// #define fast_atan __builtin_atan

static __inline float fast_sqrt(float val)  {
        union {
            int32_t tmp;
            float val;
        } u;
    
        u.val = val;
        // Remove last bit so 1.0 gives 1.0
        u.tmp -= 1<<23;
        // tmp is now an approximation to logbase2(val)
        u.tmp >>= 1; // divide by 2
        u.tmp += 1<<29; // add 64 to exponent: (e+127)/2 =(e/2)+63,
        // that represents (e/2)-64 but we want e/2
        return u.val;
}


static __inline float fast_sqrt_2(const float x)
{
    const float xhalf = 0.5f*x;
    union {
        float x;
        int32_t i;
    } u;
    u.x = x;
    u.i = 0x5f3759df - (u.i >> 1);  // gives initial guess y0
    return x*u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
}



static __inline int F_ABS(int v) { return std::abs(v); }
static __inline double F_ABS(double v) { return std::fabs(v); }
static __inline float F_ABS(float v) { return std::fabs(v); }
static __inline auto F_ABS(auto v) { return simd::fabs(v); }

static __inline mssFloat2 F_COPYSIGN(mssFloat2 mag, mssFloat2 sign) { return simd::copysign(mag, sign); }
static __inline mssFloat8 F_COPYSIGN(float mag, mssFloat8 sign) { return simd::copysign((mssFloat8)const1(mag), sign); }
static __inline mssFloat4 F_COPYSIGN(float mag, mssFloat4 sign) { return simd::copysign((mssFloat4)const1(mag), sign); }
static __inline float4v F_COPYSIGN(float4v mag, float4v sign) { return simd::copysign(mag, sign); }
static __inline float8v F_COPYSIGN(float8v mag, float8v sign) { return simd::copysign(mag, sign); }
static __inline double4v F_COPYSIGN(double4v mag, double4v sign) { return simd::copysign(mag, sign); }
static __inline double8v F_COPYSIGN(double8v mag, double8v sign) { return simd::copysign(mag, sign); }

static __inline double8v F_MIN(double8v a, double8v b) {return simd::min(a,b); }
static __inline double4v F_MIN(double4v a, double4v b) { return simd::min(a,b); }
static __inline double2v F_MIN(double2v a, double2v b) { return simd::min(a,b); }
static __inline float8v F_MIN(float8v a, float8v b) { return simd::fmin(a,b); }
static __inline float4v F_MIN(float4v a, float4v b) { return simd::fmin(a,b); }
static __inline float2v F_MIN(float2v a, float2v b) { return simd::fmin(a,b); }
static __inline double F_MIN(float a, double b) { return std::min((double)a, b); }
static __inline float F_MIN(float a, float b) { return std::min(a, b); }
static __inline int F_MIN(int a, int b) { return std::min(a, b); }

static __inline int8v F_MAX(int8v a, int8v b) {return simd::max(a,b); }
static __inline int4v F_MAX(int4v a, int4v b) {return simd::max(a,b); }
static __inline int2v F_MAX(int2v a, int2v b) {return simd::max(a,b); }
static __inline double8v F_MAX(double8v a, double8v b) {return simd::max(a,b); }
static __inline double4v F_MAX(double4v a, double4v b) { return simd::max(a,b); }
static __inline double2v F_MAX(double2v a, double2v b) { return simd::max(a,b); }
static __inline float8v F_MAX(float8v a, float8v b) {return simd::max(a,b); }
static __inline float4v F_MAX(float4v a, float4v b) { return simd::max(a,b); }
static __inline float2v F_MAX(float2v a, float2v b) { return simd::max(a,b); }
static __inline double F_MAX(double a, double b) { return std::max(a, b); }
static __inline double F_MAX(double a, float b) { return std::max(a, (double)b); }
static __inline double F_MAX(float a, double b) { return std::max((double)a, b); }
static __inline float F_MAX(float a, float b) { return std::max(a, b); }
static __inline int F_MAX(int a, int b) { return std::max(a, b); }

static __inline float8v F_CLAMP(float8v x, float8v a, float8v b) { return simd::clamp(x, a, b); }
static __inline float4v F_CLAMP(float4v x, float4v a, float4v b) { return simd::clamp(x, a, b); }
static __inline float2v F_CLAMP(float2v x, float2v a, float2v b) { return simd::clamp(x, a, b); }
static __inline double8v F_CLAMP(double8v x, double8v a, double8v b) { return simd::clamp(x, a, b); }
static __inline double4v F_CLAMP(double4v x, double4v a, double4v b) { return simd::clamp(x, a, b); }
static __inline double2v F_CLAMP(double2v x, double2v a, double2v b) { return simd::clamp(x, a, b); }
static __inline double F_CLAMP(double x, double a, double b) { return std::clamp(x, a, b); }
static __inline double F_CLAMP(double x, float a, float b) { return std::clamp(x, (double)a, (double)b); }
static __inline double F_CLAMP(double x, float a, double b) { return std::clamp(x, (double)a, b); }
static __inline double F_CLAMP(double x, double a, float b) { return std::clamp(x, a, (double)b); }
static __inline float F_CLAMP(float x, float a, float b) { return std::clamp(x, a, b); }
static __inline float F_CLAMP(float x, double a, float b) { return std::clamp(x, (float)a, b); }
static __inline float F_CLAMP(float x, float a, double b) { return std::clamp(x, a, (float)b); }

static __inline auto F_SIN(auto s) { return simd::sin(s); }
static __inline double F_SIN(double s) { return std::sin(s); }
static __inline float F_SIN(float s) { return std::sin(s); }

static __inline int __float_as_int(float in) {
     union fi { int i; float f; } conv;
     conv.f = in;
     return conv.i;
}

static __inline float __int_as_float(int a)
{
    union {int a; float b;} u;
    u.a = a;
    return u.b;
}

__inline static double fast_logf(double a)
{
    double m, r, s, t, i, f;
    int32_t e;

 //   return __builtin_log(a);
    
    if ((a > 0.0) && (a <= 3.40e+38)) { // 0x1.fffffep+127
        m = frexpf(a, &e);
        if (m < 0.666666667) {
            m = m + m;
            e = e - 1;
        }
        i = (float)e;
        /* m in [2/3, 4/3] */
        f = m - 1.0f;
        s = f * f;
        /* Compute log1p(f) for f in [-1/3, 1/3] */
        r = fmaf(-0.130187988, f, 0.140889585); // -0x1.0aa000p-3, 0x1.208ab8p-3
        t = fmaf(-0.121489584, f, 0.139809534); // -0x1.f19f10p-4, 0x1.1e5476p-3
        r = fmaf(r, s, t);
        r = fmaf(r, f, -0.166845024); // -0x1.55b2d8p-3
        r = fmaf(r, f, 0.200121149); //  0x1.99d91ep-3
        r = fmaf(r, f, -0.249996364); // -0x1.fffe18p-3
        r = fmaf(r, f, 0.333331943); //  0x1.5554f8p-2
        r = fmaf(r, f, -0.500000000); // -0x1.000000p-1
        r = fmaf(r, s, f);
        r = fmaf(i, 0.693147182, r); //   0x1.62e430p-1 // log(2)
        return r;
    }
    return 0.0;
}

__inline static double fast_expf (double a)
{
  //  return __builtin_exp(a);
    
    double f, r, j, s, t;
    long i, ia;

    // exp(a) = 2**i * exp(f); i = rintf (a / log(2))
    j = fmaf (1.442695, a, 12582912.) - 12582912.; // 0x1.715476p0, 0x1.8p23
    f = fmaf (j, -6.93145752e-1, a); // -0x1.62e400p-1  // log_2_hi
    f = fmaf (j, -1.42860677e-6, f); // -0x1.7f7d1cp-20 // log_2_lo
    i = (int)j;
    // approximate r = exp(f) on interval [-log(2)/2, +log(2)/2]
    r =             1.37805939e-3;  // 0x1.694000p-10
    r = fmaf (r, f, 8.37312452e-3); // 0x1.125edcp-7
    r = fmaf (r, f, 4.16695364e-2); // 0x1.555b5ap-5
    r = fmaf (r, f, 1.66664720e-1); // 0x1.555450p-3
    r = fmaf (r, f, 4.99999851e-1); // 0x1.fffff6p-2
    r = fmaf (r, f, 1.00000000e+0); // 0x1.000000p+0
    r = fmaf (r, f, 1.00000000e+0); // 0x1.000000p+0
    // exp(a) = 2**i * r
    ia = (i > 0) ?  0 : 0x83000000;
    s = __int_as_float (0x7f000000 + ia);
    t = __int_as_float ((i << 23) - ia);
    r = r * s;
    r = r * t;
    // handle special cases: severe overflow / underflow
    if (F_ABS (a) >= 104.0) r = s * s;
    return r;
}
static __inline double F_LOG(double g) { return fast_logf(g); }
static __inline float F_LOG(float g) { return fast_logf(g); }


static __inline double fast_powf(double a, double b) {
  // return powf( a, b );
   return fast_expf(b * F_LOG(a));
}
static __inline double db2lin(double db){ // dB to linear
   return fast_powf(10.0, 0.05 * db);
}

static __inline double fast_atan( double x )
{
    return __builtin_atan( x );
    
    double a, z, p, r, q, s, t;
    // argument reduction:
    //   arctan (-x) = -arctan(x);
    //   arctan (1/x) = 1/2 * pi - arctan (x), when x > 0
    
    z = F_ABS (x);
    a = (z > 1.0) ? (1.0 / z) : z;
    s = a * a;
    q = s * s;
    // core approximation: approximate atan(x) on [0,1]
    p =            -2.0258553044340116e-5;  // -0x1.53e1d2a258e3ap-16
    t =             2.2302240345710764e-4;  //  0x1.d3b63dbb6167ap-13
    p = fma (p, q, -1.1640717779912220e-3); // -0x1.312788ddde71dp-10
    t = fma (t, q,  3.8559749383656407e-3); //  0x1.f9690c824aaf1p-9
    p = fma (p, q, -9.1845592187222193e-3); // -0x1.2cf5aabc7dbd2p-7
    t = fma (t, q,  1.6978035834594660e-2); //  0x1.162b0b2a3bcdcp-6
    p = fma (p, q, -2.5826796814492296e-2); // -0x1.a7256feb6f841p-6
    t = fma (t, q,  3.4067811082715810e-2); //  0x1.171560ce4a4ecp-5
    p = fma (p, q, -4.0926382420509999e-2); // -0x1.4f44d841450e8p-5
    t = fma (t, q,  4.6739496199158334e-2); //  0x1.7ee3d3f36bbc6p-5
    p = fma (p, q, -5.2392330054601366e-2); // -0x1.ad32ae04a9fd8p-5
    t = fma (t, q,  5.8773077721790683e-2); //  0x1.e17813d669537p-5
    p = fma (p, q, -6.6658603633512892e-2); // -0x1.11089ca9a5be4p-4
    t = fma (t, q,  7.6922129305867892e-2); //  0x1.3b12b2db5173cp-4
    p = fma (p, s, t);
    p = fma (p, s, -9.0909012354005267e-2); // -0x1.745d022f8dc5fp-4
    p = fma (p, s,  1.1111110678749421e-1); //  0x1.c71c709dfe925p-4
    p = fma (p, s, -1.4285714271334810e-1); // -0x1.2492491fa1742p-3
    p = fma (p, s,  1.9999999999755005e-1); //  0x1.99999999840cdp-3
    p = fma (p, s, -3.3333333333331838e-1); // -0x1.5555555555448p-2
    p = fma (p * s, a, a);
    // back substitution in accordance with argument reduction //
    // double-precision factorization of PI/2 courtesy of Tor Myklebust //
    r = (z > 1.0) ? fma (0.93282184640716537, 1.6839188885261840, -p) : p;
    return copysign (r, x);
}

static __inline double F_SQRT(double s) { return std::sqrt(s); }
static __inline float F_SQRT(float s) { return std::sqrt(s); }
static __inline auto F_SQRT(auto s) { return simd::sqrt(s); }

static __inline double F_ATAN(double s) { return std::atan(s); }
static __inline float F_ATAN(float s) { return std::atan(s); }
static __inline auto F_ATAN(auto s) { return simd::atan(s); }

static __inline double F_ACOS(double s) { return acos(s); }
static __inline auto F_ACOS(auto s) { return simd::acos(s); }

static __inline size_t F_MAX(size_t a, size_t b) { return std::max(a, b); }
static __inline size_t F_MIN(size_t a, size_t b) { return std::min(a, b); }
static __inline int F_CLAMP(int x, int a, int b) { return std::clamp(x, a, b); }

static __inline float __builtin_reduce_avg(const float8v &a) {
    float x = simd_reduce_add(a) * 0.5f;
    return (x != x) ? 0.0f : x;
}
static __inline float __builtin_reduce_avg(const float4v &a) {
    float x = simd_reduce_add(a) * 0.25f;
    return (x != x) ? 0.0 : x;
}
static __inline float __builtin_reduce_avg(const float2v &a) {
    float x = simd_reduce_add(a) * 0.125;
    return (x != x) ? 0.0 : x;
}

static __inline double __builtin_reduce_avg(const double8v &a) {
    double x = simd_reduce_add(a) * 0.5;
    return (x != x) ? 0.0 : x;
}
static __inline double __builtin_reduce_avg(const double4v &a) {
    double x = simd_reduce_add(a) * 0.25;
    return (x != x) ? 0.0 : x;
}
static __inline double __builtin_reduce_avg(const double2v &a) {
    double x = simd_reduce_add(a) * 0.125;
    return (x != x) ? 0.0 : x;
}

#define F_SQR(x) ((x) * (x))

static __inline double fast_isqrt(double v)
{
    float x=(float)v,xhalf = 0.5f*x;
    int32_t i = *(int32_t*)&x;
    i = 0x5f3759df-(i>>1);
    x = *(float*)&i;
    x = x*(1.5f-(xhalf*x*x));
    return (double)x;
}

static constexpr double G1 = 0.232829;
static constexpr double G2 = 0.612779;
static constexpr double G3 = 0.843573;
static constexpr double G4 = 0.942636;
static constexpr double G5 = 0.980351;
static constexpr double G6 = 0.995302;

template <typename T>
static __inline T PHASESHIFT0(const T &in, T * co)
{
    co[0]  = (co[6]  - in)    * G1 + co[4];
    co[1]  = (co[8]  - co[0]) * G3 + co[6];
    co[2]  = (co[10] - co[1]) * G5 + co[8];
    co[10] = co[9];
    co[9]  = co[2];
    co[8]  = co[7];
    co[7]  = co[1];
    co[6]  = co[5];
    co[5]  = co[0];
    co[4]  = co[3];
    co[3]  = in;
    return co[2];
}

template <typename T>
static T PHASESHIFT90(const T &in, T * co)
{
    T out;
    co[0]  = (co[6]  - in)    * G1 + co[4];
    co[1]  = (co[8]  - co[0]) * G3 + co[6];
    co[2]  = (co[10] - co[1]) * G5 + co[8];
    out = co[9];
    co[10] = co[9];
    co[9]  = co[2];
    co[8]  = co[7];
    co[7]  = co[1];
    co[6]  = co[5];
    co[5]  = co[0];
    co[4]  = co[3];
    co[3]  = in;
    return out;
}

class DenormalsOff {
   
public:
   DenormalsOff() {
       // Store current FPCR
       asm volatile("mrs %0, fpcr" : "=r" (fpcr_));
       // Set FZ (Flush-to-zero mode) bit to disable denormals
       asm volatile("orr %0, %0, (1 << 24)" : "+r" (fpcr_));
       // Write back the modified FPCR
       asm volatile("msr fpcr, %0" :: "r" (fpcr_));
   }

   ~DenormalsOff() {
       // Restore the original FPCR
       asm volatile("msr fpcr, %0" :: "r" (fpcr_));
   }
private:
   uint64_t fpcr_;
};

static __inline double fast_log2(double x) { return std::log2(x); }
static __inline double fast_log(double x) { return std::log(x); }

static __inline double F_POW(double a, double b) { return fast_powf(a,b); }
static __inline float F_POW(float a, float b) { return fast_powf(a, b); }
static __inline double8v F_POW(double8v a, double8v b) { return simd::pow(a, b); }
static __inline double4v F_POW(double4v a, double4v b) { return simd::pow(a, b); }
static __inline double2v F_POW(double2v a, double2v b) { return simd::pow(a, b); }

static __inline auto F_EXP(auto g) { return simd::exp(g); }
static __inline double F_EXP(double g) { return std::exp(g); }
static __inline float F_EXP(float g) { return std::expf(g); }

static __inline auto F_ATAN2(auto a, auto b) { return simd::atan2(a, b); }
static __inline zfloat F_ATAN2(zfloat a, zfloat b) { return __builtin_atan2(a, b); }

static __inline auto F_COS(auto g) { return simd::cos(g); }
static __inline double F_COS(double g) { return std::cos(g); }
static __inline float F_COS(float g) { return std::cos(g); }

static __inline double F_LOG10(double g) { return __builtin_log10(g); }
static __inline auto F_LOG10(auto g) { return simd::log10(g); }

static __inline mssFloat2 from_dB(mssFloat2 x) { return F_EXP( x * 0.1151292546497023 ); }
static __inline mssFloat8 from_dB(mssFloat8 x) { return F_EXP( x * 0.1151292546497023 ); }
static __inline double from_dB(double x) { return F_EXP( x * 0.1151292546497023 ); }

static __inline auto to_dB(auto g) {return (20.0f * F_LOG10(g)); }
static __inline float to_dB(float g) { return (20.0f * F_LOG10(g)); }

#define USE_NO_FINITE_MATH // -fno-finite-math-only

#if defined(USE_NO_FINITE_MATH)

static __inline bool ISNAN(float x)  { return std::isnan(x); } //return x != x; }
static __inline bool ISNAN(double x) { return std::isnan(x); } //return x != x; }

static __inline auto ISNAN(auto x) { return simd::isnan(x); }

static __inline bool ISINF(double x) { return std::isinf(x); }
static __inline bool ISINF(float x)  { return std::isinf(x); }

static __inline auto ISINF(auto x) { return simd::isnan(x); }

__inline double F_RINT(double d) {
    d += 6755399441055744.0;
    return (double)reinterpret_cast<int&>(d);
}

static __inline auto F_RINT(auto x) { return simd::rint(x); }
/*static __inline float4v F_RINT(float4v x) { return simd::rint(x); }
static __inline float2v F_RINT(float2v x) { return simd::rint(x); }

static __inline double8v F_RINT(double8v x) { return simd::rint(x); }
static __inline double4v F_RINT(double4v x) { return simd::rint(x); }
static __inline double2v F_RINT(double2v x) { return simd::rint(x); }
*/

static __inline bool ISNORM(double v) { return !ISINF(v) && !ISNAN(v) && F_ABS(v) >= 1.175494351E-38; }

#else

static __inline bool ISINF(float x) { return std::isinf(x); }
static __inline bool ISINF(double x) { return std::isinf(x); }

static __inline bool ISNAN(float x) { return std::isnan(x); }
static __inline bool ISNAN(double x) { return std::isnan(x); }

static __inline int ISNORMf(float v) { return std::isnormal(v); }
static __inline int ISNORMd(double v) { return std::isnormal(v); }

#endif



__inline double sanitize_denormal(const double &val)
{
    return (ISINF(val) || ISNAN(val)) ? 0.0 : val;
}

__inline bool is_normal(const double &val)
{
    return !ISINF(val) && !ISNAN(val);
}

#define IS_NORMAL_DECL(TYPE) \
__inline bool is_normal(const TYPE &val) \
{ \
    TYPE isnan = convertvector( (-ISNAN(val)) + (-ISINF(val)), TYPE ); \
    return (__builtin_reduce_max( isnan ) == 0.0); \
}

IS_NORMAL_DECL(float8v)
IS_NORMAL_DECL(float4v)
IS_NORMAL_DECL(float2v)
IS_NORMAL_DECL(double8v)
IS_NORMAL_DECL(double4v)
IS_NORMAL_DECL(double2v)

#define SANITIZE_DENORMAL_DECL(TYPE) \
__inline TYPE sanitize_denormal(const TYPE &val) \
{ \
    return (is_normal( val )) ? val : 0.0; \
}

SANITIZE_DENORMAL_DECL(float8v)
SANITIZE_DENORMAL_DECL(float4v)
SANITIZE_DENORMAL_DECL(float2v)
SANITIZE_DENORMAL_DECL(double8v)
SANITIZE_DENORMAL_DECL(double4v)
SANITIZE_DENORMAL_DECL(double2v)

template <typename T>
static __inline int sanitize_denormals( T *out, int count, const char *str = "" )
{
    int present = 0;
    
    while (count--) {
        T s = *out;
        if (is_normal(s))
            *out++ = s;
        else {
            *out++ = 0.0;
            present++;
        }
    }
 //   if (present > 0)
   //     LOGGER("%s: baddies found: %d", str, present);
    return present;
}


class PeakDistAnalyzer
{
    static const int MAX_DIST = 30;
public:
    PeakDistAnalyzer() {
        peak_idx = idx = 0;
    }
    
    void reset() {

        for (int i=0; i<peak_idx; i++)
            peaks[i].peak = 0.0;
    }
    
    
    void print(void) {
        
        for (int i=0; i<peak_idx; i++)
            fprintf(stderr, "[%2d] %5s : %f\n", i+1, peaks[i].name, (float)peaks[i].peak);
    }
    
    template <typename T>
    void calcGain( T *out, int count, const char *str)
    {
        double peak = 0.0;
        
        if (!str) return;
        
        if (!strncasecmp(str, peaks[0].name, 255))
            idx = 0;
        
        
        while (count--) {
            T s = *out++;
            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
                peak = F_MAX(F_ABS(s), peak);
            else
                peak = F_MAX(__builtin_reduce_max(F_ABS(s)), peak);
        }
       
        peak_idx = F_MAX(peak_idx, idx + 1);
        peaks[idx].peak = F_MAX(peaks[idx].peak, peak);
        strncpy(peaks[idx].name, str, 255);
        
        idx++;
    }

protected:
    
    struct {
        double peak;
        char name[255];
    } peaks[MAX_DIST];
    
    int idx, peak_idx;
};
#define P_STR(x) #x
#define P_STRINGIFY(x) P_STR(x)

#define PRINT_D1_CMPLX(v) fprintf(stderr, "%s(D1): [%e:%e]\n", P_STRINGIFY(v), v.re, v.im);
#define PRINT_D1(v) fprintf(stderr, "%s(D1): [%e]\n", P_STRINGIFY(v), v);
#define PRINT_D2(v) fprintf(stderr, "%s(D2): [%.3f %.3f]\n", P_STRINGIFY(v), v[0], v[1]);
#define PRINT_D8(v) fprintf(stderr, "%s(D8): [%e %e %e %e %e %e %e %e]\n", P_STRINGIFY(v), v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
#define PRINT_D4(v) fprintf(stderr, "%s(D4): [%e %e %e %e]\n", P_STRINGIFY(v), v[0], v[1], v[2], v[3]);
#define PRINT_I4(v) fprintf(stderr, "%s(I4): [%d %d %d %d]\n", P_STRINGIFY(v), v[0], v[1], v[2], v[3]);
#define PRINT_I8(v) fprintf(stderr, "%s(I8): [%d %d %d %d %d %d %d %d]\n", P_STRINGIFY(v), v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);

#define getName(var)            #var
#define CHECK_NULL(x)           if (!(x)) { fprintf(stderr, "CHECK_NULL: %s\n",getName(x)); raise(SIGINT); }

static constexpr float SQ2_2 = 0.7071067811865475244;


#define CMPLX_MAG(X,Y)          F_SQRT( (X) * (X) + (Y) * (Y) )

#ifndef D_CMPLX_T
#define D_CMPLX_T

template <typename T>
struct cmplxT {
public:
    cmplxT(T r, T i) { re = r; im = i; }
    cmplxT(const cmplxT& v) {
        re = v.re;
        im = v.im;
    }
    cmplxT(long double v) {
        re = v; im = v;
    }
    cmplxT() { re = 0.0; im = 0.0; }
    
    T mag() const {
        return F_SQRT(re * re + im * im);
    }
    
    inline cmplxT<T> operator * (const double d) {
        return cmplxT(re * d, im * d);
    }
    
    inline cmplxT<T> operator * (const cmplxT<T>& d) {
        return cmplxT(re * d.re, im * d.im);
    }
    
    inline cmplxT<T> operator + (const cmplxT<T>& d) {
        return cmplxT(re + d.re, im + d.im);
    }
    
    inline cmplxT<T> operator - (const cmplxT<T>& d) {
        return cmplxT(re - d.re, im - d.im);
    }
    
    inline void operator *= (long double d) {
        re *= d;
        im *= d;
    }
    
    inline void operator *= (const cmplxT<T> &x)
    {
        *this = cmplxT<T>(re * x.re - im * x.im, re * x.im + im * x.re);
    }
    
    inline cmplxT<T> operator / (const cmplxT<T> &x)
    {
        T denum = (x.re * x.re + x.im * x.im);
        return cmplxT<T>((re * x.re + im * x.im) / denum, (im * x.re - re * x.im) / denum);
    }
    
    inline cmplxT<T> operator ^ (const cmplxT<T> &d) {
        return cmplxT<T>(re * d.re - im * d.im, re * d.im + im * d.re);
    }

    inline cmplxT<T> operator & (const int8v& mask) {
        return cmplxT<T>(shufflevector(re,mask), shufflevector(im,mask));
    }

    inline void operator ^= (const cmplxT<T> &d) {
       *this = cmplxT<T>(re * d.re - im * d.im, re * d.im + im * d.re);
    }
    
    inline void operator += (const cmplxT<zfloat4> &d) {
        re += d.re;
        im += d.im;
    }
    
    inline void operator += (const cmplxT<zfloat8> &d) {
        re += d.re;
        im += d.im;
    }
    
    inline void operator += (const cmplxT<zfloat> &d) {
        re += d.re;
        im += d.im;
    }

    inline T real() const { return re; }
    inline T imag() const { return im; }

    T re;
    T im;
    
}; // __attribute__((packed));

#endif //D_CMPLX_T

/* SSE and co like 16-bytes aligned pointers
 * with a 64-byte alignment, we are even aligned on L2 cache lines... */

#define MALLOC_V4SF_ALIGNMENT   64
#define VALIGNED_ID             'VLGN'

static void *Valigned_malloc(size_t nb_bytes)
{
    void *p, *p0 = malloc(nb_bytes + MALLOC_V4SF_ALIGNMENT);
    if (!p0) return (void *) 0;
    p = (void *) (((size_t) p0 + MALLOC_V4SF_ALIGNMENT) & (~((size_t) (MALLOC_V4SF_ALIGNMENT-1))));
    *(uint32_t*)((void **) p - 1) = VALIGNED_ID;
    *((void **) p - 2) = p0;
    return p;
}

static void Valigned_free(void *p)
{
    if (p && *(uint32_t*)((void **) p - 1) == VALIGNED_ID) {
        void *ptr = *((void**) p - 2);
      //  int64_t sz = bigpool.malloc_size(ptr);
        free( ptr );
      //  if ((sz/1048576L) == 0)
      //      fprintf(stderr, "Valigned_free( %3lld kb)\n", (int64_t)(sz/1024));
      //  else
      //      fprintf(stderr, "Valigned_free( %3lld MB )\n", (int64_t)(sz/1048576L));
    }
}

static inline void *Valigned_ptr(void *p)
{
    return (p) ? *((void **) p - 1) : NULL;
}

template <typename T>
static inline T *callocT( int64_t k )
{
    T *p = (T*)Valigned_malloc (k * sizeof(T));
    if (!p)
        return NULL;
    memset (p, 0, k * sizeof (T));
    return p;
}

static inline void callocT_free( void *p )
{
    Valigned_free(p);
}

#define THROW_MEM(cond)  if (cond) throw (Converror (Converror::MEM_ALLOC));
