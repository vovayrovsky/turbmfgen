/**

    \file
    \brief rk4 lib for test particle

    solve equation y' = f(t, y), where t is time and y - is {rx, ry, rz, vx, vy, vz} 

*/

#ifndef RK4_GUARD
#define RK4_GUARD

#include <math.h>

#include "v3d.h"
#include "SI_units.h"

/**

    \brief variables structure for rk4

*/
typedef struct
    {
    v3d r;
    v3d v;
    double Energy;
    } var_list;


/**

    \brief External function that should return magnetic field value at the position

    \param[in] radius vector

    \return magnetic field at this point
*/
extern v3d GetField (v3d);


/**

    \brief Get derrivatives

    \param[in] t time 
    \param[in] y variables list 
    \param[in] mass mass of particle
    \param[in] charge charge of particles

    \return magnetic field at this point
*/
var_list derivatives (double t, var_list y, double mass, double charge);

/**

    \brief rk4 step

    \param[in] prev previous state 
    \param[in] mass mass of particle
    \param[in] charge of particles
    \param[in] dt delta t
    \param[in] t  current time

    \return new state
*/
var_list step (var_list prev, double mass, double charge, double dt, double t);

#endif
