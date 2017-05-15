#pragma once
#ifndef TestUtilities_hpp__c71466562e3e42a28053db4f68d636b9
#define TestUtilities_hpp__c71466562e3e42a28053db4f68d636b9

#include <catch.hpp>
#include <iomanip>
#include <iostream>

#include <colorsystem.hpp>

struct putFloat_
{
    const float v_;
    explicit putFloat_(float v) : v_{v}
    {
        /* NO-OP */
    }
    std::ostream &write(std::ostream &o) const
    {
        auto f = o.setf(std::ios::fixed, std::ios::floatfield);
        o.width(12);
        auto prec = o.precision(6);
        o << v_;
        o.precision(prec);
        o.setf(f, std::ios::floatfield);
        return o;
    }
};

inline std::ostream &operator<<(std::ostream &o, const putFloat_ &value)
{
    return value.write(o);
}

inline putFloat_ putFloat(float v)
{
    return putFloat_{v};
}

inline std::ostream &writeTo(std::ostream &o, const ColorSystem::Vector3 &v)
{
    o << "(" << putFloat(v[0]) << " " << putFloat(v[1]) << " " << putFloat(v[2]) << ")";
    return o;
}

inline std::ostream &operator<<(std::ostream &o, const ColorSystem::Vector3 &v)
{
    return writeTo(o, v);
}

inline std::ostream &writeTo(std::ostream &o, const ColorSystem::Matrix3 &m)
{
    o << "|" << putFloat(m[0]) << " " << putFloat(m[1]) << " " << putFloat(m[2]) << "|\n"
      << "|" << putFloat(m[3]) << " " << putFloat(m[4]) << " " << putFloat(m[5]) << "|\n"
      << "|" << putFloat(m[6]) << " " << putFloat(m[7]) << " " << putFloat(m[8]) << "|\n";
    return o;
}

inline std::ostream &operator<<(std::ostream &o, const ColorSystem::Matrix3 &m)
{
    return writeTo(o, m);
}

inline std::ostream &writeTo(std::ostream &o, const ColorSystem::Tristimulus &s)
{
    return writeTo(o, s.v_);
}

inline std::ostream &operator<<(std::ostream &o, const ColorSystem::Tristimulus &s)
{
    return writeTo(o, s);
}

namespace Catch
{
/// ColorSpace::Vector3 specialized string converter.
template <>
struct StringMaker<ColorSystem::Vector3>
{
    static std::string convert(const ColorSystem::Vector3 &v)
    {
        std::ostringstream O;
        writeTo(O, v);
        return O.str();
    }
};
/// ColorSpace::Matrix3 specialized string converter.
template <>
struct StringMaker<ColorSystem::Matrix3>
{
    static std::string convert(const ColorSystem::Matrix3 &v)
    {
        std::ostringstream O;
        writeTo(O, v);
        return O.str();
    }
};

/// ColorSpace::Tristimulus specialized string converter.
template <>
struct StringMaker<ColorSystem::Tristimulus>
{
    static std::string convert(const ColorSystem::Tristimulus &v)
    {
        std::ostringstream O;
        writeTo(O, v);
        return O.str();
    }
};
}

template <typename T_>
class ApproxEquals final : public Catch::MatcherBase<T_>
{
    const double epsilon_ = std::numeric_limits<float>::epsilon() * 100;
    const T_ &   expected_;

  public:
    ApproxEquals(const T_ &expected, double eps)
        : expected_{expected}, epsilon_{eps}
    {
        /* NO-OP */
    }

    bool match(const T_ &value) const override
    {
        if (value.size() != expected_.size())
        {
            return false;
        }
        for (int_fast32_t i = 0; i < value.size(); ++i)
        {
            if (value[i] != Approx(expected_[i]).epsilon(epsilon_))
            {
                return false;
            }
        }
        return true;
    }

    std::string describe() const override
    {
        std::ostringstream O;
        O << "matches to\n"
          << expected_ << " with epsilon = " << epsilon_;
        return O.str();
    }
};

inline ApproxEquals<ColorSystem::Matrix3> IsApproxEquals(const ColorSystem::Matrix3 &m, double eps)
{
    return {m, eps};
}

inline ApproxEquals<ColorSystem::Vector3> IsApproxEquals(const ColorSystem::Vector3 &v, double eps)
{
    return {v, eps};
}

inline ApproxEquals<ColorSystem::Tristimulus> IsApproxEquals(const ColorSystem::Tristimulus &s, double eps)
{
    return {s, eps};
}
#endif /* TestUtilities_hpp__c71466562e3e42a28053db4f68d636b9 */
