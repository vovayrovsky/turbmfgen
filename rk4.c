#include "rk4.h"
#include "stdio.h"

var_list derivatives (double t, var_list y, double mass, double charge)
    {
    v3d velocity = mul (normalize (y.v), light_c);    
    v3d Bfield = GetField (y.r);

    v3d dv_over_dt = mul (vec_mul (velocity, Bfield), charge * light_c * light_c / y.Energy); 
    var_list retval = {velocity, dv_over_dt};

    return retval;
    }  
    
const double CashKarpC [6] = {0., 0.2, 0.3, 0.6, 1., 0.875};
//const double CashKarpB [6] = {2825./27648., 0., 18575./48384., 13525./55296., 277./.14336, 1./4.};
const double CashKarpB [6] = {37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.};

const double CashKarpA[] = {
	0., 0., 0., 0., 0., 0.,
	1. / 5., 0., 0., 0., 0., 0.,
	3. / 40., 9. / 40., 0., 0., 0., 0.,
	3. / 10., -9. / 10., 6. / 5., 0., 0., 0.,
	-11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0.,
	1631. / 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253. / 4096., 0.
    };

var_list step (var_list prev, double mass, double charge, double dt, double t) 
    {
    var_list k[6] = {};

    for (int i = 0; i < 6; i++)
        {
        var_list delta = {prev.r, prev.v, prev.Energy};
        for (int j = 0; j < i; j++)
            {
            delta.r = add (delta.r, mul (k[j].r, dt * CashKarpA[i*6 + j]));
            delta.v = add (delta.v, mul (k[j].v, dt * CashKarpA[i*6 + j]));
            }

        k[i] = derivatives (t, delta, mass, charge);
        }
    
    var_list retval = {prev.r, prev.v, prev.Energy};
    for (int i = 0; i < 6; i++)
        {
        retval.r = add (retval.r, mul (k[i].r, CashKarpB[i] * dt));
        retval.v = add (retval.v, mul (k[i].v, CashKarpB[i] * dt));
        }
   
    retval.v = mul (normalize (retval.v), light_c);

    return retval;
    }

