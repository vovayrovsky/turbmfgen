/**

    \file
    \brief 3d double precision real vector lib

*/

#ifndef V3D_GUARD
#define V3D_GUARD

/**

    \brief 3d vector structure

*/
typedef struct {
    ///x component
    double x;
    ///y component
    double y;
    ///z component
    double z;
    } v3d;

/**

    \brief adding two vectors

    \param[in] a first vector
    \param[in] b second vector

    \return sum of a and b

*/
v3d add (v3d a, v3d b);

/**

    \brief substracting two vectors

    \param[in] a first vector
    \param[in] b second vector

    \return difference between a and b

*/
v3d sub (v3d a, v3d b);


/**

    \brief multiplying vector by number

    \param[in] a vector
    \param[in] b number

    \return b * a

*/
v3d mul (v3d a, double b);

/**

    \brief dividing vector by number

    \param[in] a vector
    \param[in] b number

    \return a / b

*/
v3d div (v3d a, double b);


/**

    \brief vector multiplication

    \param[in] a first vector
    \param[in] b second vector

    \return vector product of a and b

*/
v3d vec_mul (v3d a, v3d b);

/**

    \brief substracting two vectors

    \param[in] a first vector
    \param[in] b second vector

    \return scalar product of a and b

*/
double scl_mul (v3d a, v3d b);

/**

    \brief length of vector

    \param[in] a vector

    \return length of a

*/
double module (v3d a);

/**

    \brief normalizing vector

    \param[in] a vector

    \return normalized a

*/
v3d normalize (v3d a);

/**

    \brief construct vector from it's polar coordinates

    \param[in] a        length of vector
    \param[in] thetta   angle between z axis and vector
    \param[in] fi       angle between x axis and projection of vector on XY plane

    \return constructed vector

*/
v3d FromPolar (double a, double thetta, double fi);

#endif
