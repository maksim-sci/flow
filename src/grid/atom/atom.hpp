#pragma once

#include <memory>

#include "type.hpp"

#include "../../sgs.hpp"

namespace grid {
    namespace atom {
        typedef std::shared_ptr<Type> pType_t;
        class Atom {
            pType_t material;
            bool fixu;
            double q;
            double u;
        public:
            inline Atom(pType_t t):material(t),q(material->Q()),u(0),fixu(false) {};

            inline double Q() const {return q;};
            inline void Q(double _q) {q=_q;};

            inline double U() const {return u;};
            inline void U(double _u) {if(fixu) return; u=_u;};

            inline double T() const {return 300*sgs::KELVIN;/*TODO*/};
            inline void T(double _t) const {/*TODO*/};

            inline const pType_t Material() const {return material;};
            inline const void Material(pType_t mat) {material = mat;q=material->Q();};

            inline bool fixU() {return fixu;};
            inline void fixU(bool b) {fixu=b;};
        };
    }
}