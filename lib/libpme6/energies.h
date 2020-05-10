#ifndef ENERGIES_H
#define ENERGIES_H

#include <ostream>

struct energies {
  double direct;
  double recip;
  double extra;
  double self;
  double masked;

  energies()
      : direct(0.0), recip(0.0), extra(0.0),
        self(0.0), masked(0.0) {}

  double total() const {
    return direct + recip + extra - self - masked;
  }

  friend std::ostream& operator<<(std::ostream& stream, energies& ener) {
    return stream << ener.total()
                  << " = direct + recip + extra - self - masked = "
                  << ener.direct << " + " 
                  << ener.recip  << " + "
                  << ener.extra  << " - " 
                  << ener.self   << " - "
                  << ener.masked;
  }
};

#endif
