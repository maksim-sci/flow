#pragma once

namespace sgs {
    constexpr double CENTIMETER = 1.;
    constexpr double GRAM = 1.;
    constexpr double SECOND = 1.;
    //custom values
    constexpr double METER = 100.*CENTIMETER;
    constexpr double ANGSTROM = 1e-10*METER;
    constexpr double MINUTE = 60.*SECOND;
    constexpr double KILOGRAM = 1000.*GRAM;
    //complex values
    constexpr double KELVIN = 1.;
    constexpr double JOULE = KILOGRAM*(METER/SECOND)*(METER/SECOND);
    constexpr double ERG = GRAM*(CENTIMETER/SECOND)*(CENTIMETER/SECOND);
    //constants
    constexpr double PLANCK  = 1.0544e-27*ERG*SECOND;
    constexpr double ELVOLT  = 1.6e-12*ERG;
    
    constexpr double VOLT  = 0.003335640484669; //TODO fix it
    
    constexpr double ELCHARGE = 4.80320427e-10;//TODO   *f64::powf(CENTIMETER, 1.5)*GRAM.powf(0.5)/SECOND;
    constexpr double ELMASS = 9.109383e-31*KILOGRAM;//TODO   *f64::powf(CENTIMETER, 1.5)*GRAM.powf(0.5)/SECOND;

    constexpr double BOLZMAN = 5.380649e-15*ERG;//TODO   *f64::powf(CENTIMETER, 1.5)*GRAM.powf(0.5)/SECOND;
}