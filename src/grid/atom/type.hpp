#pragma once
#include <string>
#include <fmt/format.h>
namespace grid {
    namespace atom {

    using std::string;
    typedef unsigned int TypeId;
        class Type{
            double baseq;
            double v;
            unsigned int id;
            string name;

        public:
            inline Type():baseq(0),id(0),name(""),v(0) {};
            inline Type(TypeId _id):baseq(0),id(_id),name(""),v(0) {};
            inline Type(double q,TypeId _id):baseq(q),id(_id),name(""),v(0) {};
            inline Type(double q,TypeId _id,string _name):baseq(q),id(_id),name(_name),v(0) {};
            inline Type(double q=0,TypeId _id=0,string _name="",double _v=0):baseq(q),id(_id),name(_name),v(_v) {};
            

            inline bool operator==(const Type& rhs) const {return id==rhs.id;};
            inline bool operator!=(const Type& rhs) const {return !(*this==rhs);};
            inline bool operator<(const Type& rhs) const {return id<rhs.id;};
            inline bool operator>(const Type& rhs) const {return id>rhs.id;};

            inline double Q() const {return baseq;};
            inline void Q(double newQ) {baseq=newQ;};
            inline void V(double _v) {v=_v;};
            inline double V() {return v;};
            inline string Name() const {return name;};

            inline TypeId Id() const {return id;};
            


        };
    }
}

template <> struct fmt::formatter<grid::atom::Type> {
  inline constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && *it != '}') throw format_error("invalid format");
    return it;
  }
  template <typename FormatContext>
  inline auto format(const grid::atom::Type& t, FormatContext& ctx) const -> decltype(ctx.out()) {

              return fmt::format_to(ctx.out(), "atom::type(q:{},id:{},name:{})", t.Q(),t.Id(),t.Name());
  }
};