/**

    \file
    \brief Simple rand.h wrapper

*/

#ifndef RANDOM_GUARD
#define RANDOM_GUARD

///path to the random device
#define RANDOM_DEV "/dev/random"

///no error return code
#define NO_ERROR                0x000
///DEVRANDOM opening error
#define DEVRANDOM_OPEN_ERROR    0x100
///DEVRANDOM reading error
#define DEVRANDOM_READ_ERROR    0x101


/**

    \brief Initializing random number generator
    
    Initializing rand.h seed from /dev/random
       

    \retval error code (see random.h)
*/
 int InitializeRandom();

/**

    \brief Get double random number in given range
    
    \param[in] max maximum value of range
    \param[in] min minimum value of range

    \return random number
*/
double GetRandd (double max, double min);

/**

    \brief Get int random number in given range
    
    \param[in] max maximum value of range
    \param[in] min minimum value of range

    \return random number
*/
int GetRandi (int max, int min);

/**

    \brief Get int random number from RAND_MIN ro RAND_MAX (defined in stdlib.h)
    
    Actually, it's just rand() wrapper

    \return random number
*/
int GetRandiMax();

#endif

//EOF
