#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "SI_units.c"

#include "random.h"
#include "v3d.h"
#include "rk4.h"

#define FIELD_K_NUMBER 5000
#define NUM_OF_PARTICLES 32

volatile    char initialized = 0;

volatile    v3d      k       [FIELD_K_NUMBER] = {};
volatile    double   phase   [FIELD_K_NUMBER] = {};
volatile    v3d      pl      [FIELD_K_NUMBER] = {};

const double max_l = 100. * pc;
const double min_l = 100. * au; 

const double Brms = 6. * uG;
double Bb_val = 0.;

const double spectral_index = -5./3.;

double GetFieldAmplitude (double Brms, double max_k, double min_k, double k);

void InitTurbulence();
v3d GetRandomFieldLinear (v3d position);
v3d GetRandomFieldExp (v3d position);

void GenToFile (v3d size, v3d origin, int steps, FILE* fout);
v3d GetField (v3d r);

extern double CashKarpA[];
extern double CashKarpB[];
extern double CashKarpC[];

int main (int argc, char** argv)
    {
    #ifdef _OPENMP
    
    printf ("Caution: The program was compiled with OpenMP and can consume all CPU resources of your PC!\n");
    printf ("Start program with %d cores\n", omp_get_num_procs());

    #endif

    int errorCode = InitializeRandom();

    if (errorCode)
        {
        printf ("Error: %x\n", errorCode);
    
        return 0;
        }
    FILE* fout = NULL;
    
    if (argc > 1)
        {
        printf ("output: %s\n", argv[1]);
        fout = fopen (argv[1], "w");
        }
    else
        fout = fopen ("test_out.txt", "w");
   
    fprintf (fout, "#id\tL\tx\ty\tz\tDxx\tDyy\tdL\tvx\tvy\tvz\n");

    InitTurbulence();

//omp_set_nested (1);

    double Energy_val = 0.;
    FILE* fenergy = fopen ("energy", "r");
    if (fenergy == NULL)
        {
	printf ("Error! No energy file\n");

	return 1;
	}

    FILE* fBb = fopen ("Bb", "r");
    if (fBb != NULL)
        fscanf (fenergy, "%lg", &Energy_val);
    
    printf ("Energy: %lg TeV\n", Energy_val);

    fscanf (fBb, "%lg", &Bb_val);
    printf ("Bb: %lg\n", Bb_val);

    #pragma omp parallel for shared(fout, CashKarpA, CashKarpB, CashKarpC)
    for (int particle = 0; particle < NUM_OF_PARTICLES; particle++)
        {
        double Energy = Energy_val * TeV;
   
	printf ("Energy: %lg TeV\n", Energy_val);

        v3d direction = FromPolar (1., GetRandd (2. * M_PI, 0.), GetRandd (2. * M_PI, 0.));
        var_list state = {{0., 0., 0.}, mul (normalize (direction), light_c), Energy};
    
        double time = 0.;
        double L = 0.;
        double Lstep = 10. * (Energy_val > 10000? 10000 : Energy_val) * au;
    
        double dt = Lstep / light_c;
   
        //printf ("R: %lg au\n", p_mass * light_c / (e * 6. * uG) / sqrt (1. - pow (v, 2.)/pow (light_c, 2.)) / au);
        //printf ("R: %lg au\n", Energy/e/(6. * uG) / au);
        printf ("R: %lg au\n", 3e11 * (Energy / (10. * GeV)) / 6. /au);
 
        double mileage = 800000. * pc;

        long long unsigned int num = mileage / Lstep; 
        // int percentage = -1;

        printf ("num %llu\n", num);

        for (long long unsigned i = 0; i < num; i++)
            {
            double Dxx = 0.;
            double Dyy = 0.;

            /*if ((i * 100/ num) > percentage)
                {
                percentage = i * 100 / num;
                printf ("done %d %%\n", percentage);
                }
		*/
            
	    //if (time != 0.)
            //    {
            //    Dxx = pow (state.r.z, 2.) / time;
            //    Dyy = (pow (state.r.x, 2.) + pow (state.r.y, 2.)) / time;
            //    }

            //v = module(state.v) * light_c / sqrt (pow (p_mass * light_c, 2.) + pow (module(state.v), 2.));
            dt = Lstep / light_c;

            var_list next = step (state, p_mass, e, dt, time);

            v3d deltaL = sub (next.r, state.r);

            //printf ("dL: %lg au, dt: %lg s, v: %lg\n", module(deltaL) / au, dt, module(state.v));
	
            #pragma omp critical (fileout)
	    {
            if (i % 40 == 0)
	        fprintf (fout, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", 
                        particle, L / au, state.r.x / au, state.r.y / au, state.r.z / au, Dxx, Dyy, module(deltaL)/pc,
                        state.v.x, state.v.y, state.v.z);
            }
            #pragma omp flush
	    if (i % 10000 == 9999)
	        fflush (fout);

            L = (i + 1) * Lstep;
            state = next;

            time += dt;
            }
        }
    /*
    int num = 10000;

    v3d mean = {};
    v3d sqmean = {};

    int percentage = -1;

    int num_of_big = 0;

    FILE* field_out = fopen ("field_out", "w");

    for (int i = 0; i < num; i++)
        {
        if ((i * 100/ num) > percentage)
            {
            printf ("done %d %%\n", percentage);
            percentage = i * 100 / num;
            }

        v3d buf = {GetRandd (100. * pc, -100. * pc), GetRandd (100. * pc, -100. * pc), GetRandd (100. * pc, - 100. * pc)};
        //v3d buf = {i * 10. *  au, 0, 0};
        v3d field = GetRandomFieldExp (buf);
       
        fprintf (field_out, "%lg %lg %lg %lg %lg %lg\n", field.x / uG, field.y / uG, field.z / uG, 
                pow (field.x / uG, 2.), pow (field.y / uG, 2), pow (field.z / uG, 2)); 
        
        if (module (field) > 1.)
            //printf ("%lg uG\n", module (field) / uG);
            num_of_big++;

        mean = add (mean, field);
        sqmean = add (sqmean, (v3d) {
                                    field.x * field.x,
                                    field.y * field.y,
                                    field.z * field.z
                                    });
        //printf ("a[0]: %lg\n", module (pl[0]));

        }
   
    mean   = div (mean,   (double) num);
    sqmean = div (sqmean, (double) num);

    fclose (field_out);
    
    printf ("num of big: %d\n", num_of_big);
    printf ("Bmean: %lg {%lg, %lg, %lg} uG\n", module (mean) / uG, mean.x / uG, mean.y / uG, mean.z / uG);
    printf ("Brms: %lg {%lg, %lg, %lg} uG\n", sqrt (sqmean.x + sqmean.y + sqmean.z) / uG,
                                              sqrt (sqmean.x) / uG, sqrt (sqmean.y) / uG, sqrt (sqmean.z) / uG);
    */
    /**/
    //v3d field = GetRandomFieldLinear (buf);
    //printf ("B = {%.7lg, %.7lg, %.7lg} uG, |B| = %.7lg uG\n", field.x / uG, field.y / uG, field.z / uG, module(field) /uG);
    /*
    FILE* fout = fopen ("testfield.mf", "w");

    v3d size = {1. * pc, 1. * pc, 1. * pc};
    v3d origin = {0., 0., 0.};
    int steps = 100;

    GenToFile (size, origin, steps, fout);

    fclose (fout);
    
    fclose (fout);
    */
    printf ("Compiled at %s %s\n", __DATE__, __TIME__);

    return 0;
    }

v3d GetField (v3d r)
    {
    return add (GetRandomFieldExp (r), (v3d) {0., 0., Bb_val * 6 * uG});
    //return GetRandomFieldExp (r);
    //return (v3d) {0., 0., 6. * uG};
    }

void GenToFile (v3d size, v3d origin, int steps, FILE* fout)
    {
    typedef struct {
        char format1;
        char format2;

        double x_size;
        double y_size;
        double z_size;
        
        int steps;

        double x_origin;
        double y_origin;
        double z_origin;
        char format3;
        } Header;
         
    typedef struct {
        double x;
        double y;
        double z;
        } coords;

    Header head = {'M', 'F', size.x, size.y, size.z, steps, origin.x, origin.y, origin.z, 'I'}; 
#define FOUT_DEF(name)\
fwrite(& (name), sizeof(name), 1, fout)

    FOUT_DEF (head.format1);
    FOUT_DEF (head.format2);
    FOUT_DEF (size.x);
    FOUT_DEF (size.y);
    FOUT_DEF (size.z);
    FOUT_DEF (steps);
    FOUT_DEF (origin.x);
    FOUT_DEF (origin.y);
    FOUT_DEF (origin.z);
    FOUT_DEF (head.format3);

#undef FOUT_DEF

    //printf ("%d %d %d %d\n", sizeof (int), sizeof (double), sizeof(char), sizeof(Header));

    for (int x = 0; x < steps; x++)
        for (int y = 0; y < steps; y++)
            for (int z = 0; z < steps; z++)
                {
                v3d delta = {size.x / (double)steps * (double) x,
                             size.y / (double)steps * (double) y,
                             size.z / (double)steps * (double) z};
    
                v3d field = GetRandomFieldLinear (add (origin, delta));
                
                //v3d field = {1., 0., 0.};

                //fwrite(&field, sizeof(field), 1, fout); 
                }
    }

void InitTurbulence ()
    {
    if (!initialized)
        {
        double sAi = 0.;

        for (int i = 0; i < FIELD_K_NUMBER; i++)
            {
            double l = min_l * pow (max_l / min_l, (double) i /FIELD_K_NUMBER);
            double dl = log (max_l / min_l) / FIELD_K_NUMBER * l; 
            
            double module_k = 2. * M_PI/l;
            double dk = module_k / l * dl;
             
            double field = GetFieldAmplitude (Brms, 2. * M_PI/min_l, 2. * M_PI /max_l, module_k);

            sAi += field * dk;
            }

        for (int i = 0; i < FIELD_K_NUMBER; i++)
            {
            double l = min_l * pow (max_l / min_l, (double) i /FIELD_K_NUMBER);
            double dl = log (max_l / min_l) / (double)FIELD_K_NUMBER * l; 
            
            double module_k = 2. * M_PI/l;
            double dk = module_k / l * dl;
            
            double tetta = acos (GetRandd (1., -1.));
            double fi = GetRandd (2. * M_PI, 0.);
            k[i] = FromPolar (module_k, tetta, fi);
            phase[i] = GetRandd (2. * M_PI, 0);
         
            double angle = GetRandd (2. * M_PI, 0);
            double amplitude = GetFieldAmplitude (Brms, 2. * M_PI/min_l, 2. * M_PI /max_l, module_k);

            pl[i] = (v3d) {-sin (fi) * cos (angle) + cos (tetta) * cos (fi) * sin (angle),
                            cos (fi) * cos (angle) + cos (tetta) * sin (fi) * sin (angle),
                            -sin (tetta) * sin (angle)};
            pl[i] = normalize(pl[i]);
            
            pl[i] = mul (pl[i], sqrt (2. * amplitude * dk / sAi));
            }

        printf ("sAi: %lg uG\n", sAi / uG);

        initialized = 1;
        printf ("Init turbulence\n");
        }
    }

v3d GetRandomFieldExp (v3d position)
    {
    double sai = 0.;
    v3d retval = {0., 0., 0.};
   
    //#pragma omp parallel for shared(pl, k, position, phase, retval)
    for (int i = 0; i < FIELD_K_NUMBER; i++)
        {
        v3d dB = mul (pl[i], cos (scl_mul (k[i], position) + phase[i]));
        //#pragma omp critical
        retval = add (retval, dB);        
        //printf ("mpli %lg\n", module (pl[i]) / uG);
        //sai += module (pl[i]);
        //if (module(pl[i]) > 0)
        //    printf ("at %d %lg\n", i, module(pl[i]));
        }
    
    //printf ("sai %lg uG\n", sai / uG);
    //printf ("retval %lg uG\n", module (retval) / uG);

    retval = mul (retval, Brms);

    return retval;
    }


v3d GetRandomFieldLinear (v3d position)
    {
    if (!initialized)
        {
        double sAi = 0.;

        for (int i = 0; i < FIELD_K_NUMBER; i++)
            {
            double dl = 1. / FIELD_K_NUMBER * (max_l - min_l); 
            double l = (double) i * dl + min_l;
            double module_k = 2. * M_PI/l;
            double dk = module_k / l * dl;
             
            double field = GetFieldAmplitude (Brms, 2. * M_PI/min_l, 2. * M_PI /max_l, module_k);

            sAi += field * dk;
            }

        FILE* fout = fopen ("xi.txt", "w");

        for (int i = 0; i < FIELD_K_NUMBER; i++)
            {
            double dl = 1. / FIELD_K_NUMBER * (max_l - min_l); 
            double l = (double) i * dl + min_l;
            double module_k = 2. * M_PI/l;
            double dk = module_k / l * dl;
             
            double tetta = acos (GetRandd (1., -1.));
            double fi = GetRandd (2. * M_PI, 0.);
            k[i] = FromPolar (module_k, tetta, fi);
            phase[i] = GetRandd (2. * M_PI, 0);
         
            /*
            v3d kp1 = {0., 0., 0.};
            
            kp1.x = -k[i].y - k[i].z;
            kp1.y = k[i].x;
            kp1.z = k[i].x;
            
            v3d kp2 = vec_mul (k[i], kp1);

            //if (scl_mul(kp1, kp2) > 0.00001)
            //    printf ("D");
            
            kp1 = normalize (kp1);
            kp2 = normalize (kp2);
            */
            double angle = GetRandd (2. * M_PI, 0);
            double amplitude = GetFieldAmplitude (Brms, 2. * M_PI/min_l, 2. * M_PI /max_l, module_k);

            //if (module (kp1) - 1. > 0.00001)
            //    printf ("A");
            //if (module (kp2) - 1. > 0.00001)
            //    printf ("B");

            //if (scl_mul(kp1, kp2) > 0.00001)
            //    printf ("K");

            //pl[i] = add (mul (kp1, cos (angle)), mul (kp2, sin (angle))); //TODO: check it
            //printf ("%lg au ^ -1\t%lg uG\n", module_k * au, amplitude / uG);
            //if (module(pl[i]) - 1. > 0.0001)
            //    printf ("%lg ", module(pl[i]));
            
            pl[i] = (v3d) {-sin (fi) * cos (angle) + cos (tetta) * cos (fi) * sin (angle),
                            cos (fi) * cos (angle) + cos (tetta) * sin (fi) * sin (angle),
                            -sin (tetta) * sin (angle)};
            pl[i] = normalize(pl[i]);
            
            //printf ("mpli %lg\n", module (pl[i]));
   
            fprintf (fout, "%lg %lg %lg\n", pl[i].x, pl[i].y, pl[i].z);

            if (scl_mul (pl[i], k[i]) > 0.0001)
                {
                printf ("Bad polarization\n");
                exit(0);
                }
        
            //printf ("%lg %lg\n", amplitude * dk, sAi);

            pl[i] = mul (pl[i], sqrt (2. * amplitude * dk / sAi));
            }

        fclose (fout);

        printf ("sAi: %lg uG\n", sAi / uG);

        initialized = 1;
        printf ("Init turbulence\n");
        }

    double sai = 0.;
    v3d retval = {0., 0., 0.};
    
    for (int i = 0; i < FIELD_K_NUMBER; i++)
        {
        v3d pli = pl[i];
        retval = add (retval, mul (pli, cos (scl_mul (k[i], position) + phase[i])));        
        //printf ("mpli %lg\n", module (pl[i]) / uG);
        sai += module (pl[i]);
        //if (module(pl[i]) > 0)
        //    printf ("at %d %lg\n", i, module(pl[i]));
        }
    
    //printf ("sai %lg uG\n", sai / uG);
    //printf ("retval %lg uG\n", module (retval) / uG);

    retval = mul (retval, Brms);

    return retval;
    }


double GetFieldAmplitude (double Bint, double max_k, double min_k, double k)
    {
    //ln (dB/dk) = -11/3 * ln(k) + ln(alpha);
    //integral from k_min to k_max dB/dk = Bint
   
    //return Bint / (max_k - min_k);

    if (k < min_k || k > max_k)
        return 0.; 

    double alpha = Bint / ((1./(spectral_index - 1.)) * (pow (max_k, (spectral_index - 1.)) - pow (min_k, (spectral_index - 1.))));

    return alpha*pow (k, (spectral_index - 2.));
    }
