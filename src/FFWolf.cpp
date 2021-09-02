#include "FFWolf.h"

/**
 * Author: Greg Schwing
 * Description: Implementation of Wolf Method as a alternative to the Ewald Summation
 * for very large systems. 
 * Citation: Wolf, D.; Keblinski, P.; Phillpot, S. R.; Eggebrecht, J. Exact
 *   method for the simulation of Coulombic systems by spherically
 *   truncated, pairwise r −1 summation. J. Chem. Phys. 1999, 110, 8254−
 *   8282.
 */

double FF_WOLF::CalcCoulomb(const double dist, const double qi_qj_Fact,
                            const uint b, const double lambda) const{
    double wolf_electrostatic;
    wolf_electrostatic = erfc(alpha * dist)/dist;
    wolf_electrostatic -= erfc(alpha * Rc)/Rc;
    wolf_electrostatic *= qi_qj_Fact;
    return wolf_electrostatic;
}