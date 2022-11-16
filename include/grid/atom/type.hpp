#pragma once
#include <string>
#include <fmt/format.h>
namespace grid {
    namespace atom {

    using std::string;
    typedef unsigned int TypeId;
        class Type{
            double baseq;
            unsigned int id;
            string name;

        public:
            inline Type():baseq(0),id(0),name("") {};
            inline Type(TypeId _id):baseq(0),id(_id),name("") {};
            inline Type(double q,TypeId _id):baseq(q),id(_id),name("") {};
            inline Type(double q,TypeId _id,string _name):baseq(q),id(_id),name(_name) {};
            

            inline bool operator==(const Type& rhs) const {return id==rhs.id;};
            inline bool operator!=(const Type& rhs) const {return !(*this==rhs);};
            inline bool operator<(const Type& rhs) const {return id<rhs.id;};
            inline bool operator>(const Type& rhs) const {return id>rhs.id;};

            inline double Q() const {return baseq;};
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