#ifndef _MITRE_PDA_H_
#define _MITRE_PDA_H_

/* MACROS */

#define log_round(x) ( (x)/floor(x) <= ceil(x)/(x) ? floor(x) : ceil(x) )
#define min(a,b) ( (a) <= (b) ? (a) : (b) )
#define max(a,b) ( (a) >= (b) ? (a) : (b) )

/* The fact that C arrays' index start from 0 is often inconvenient in the context of
 * this library, making it necessary to subtract 1 from the iterator variable before
 * accessing the array position.
 * By using the slack variable it is possible to avoid this and possibly obtain some
 * performance gain, although wasting an irrelevant O(1) amount of memory.
 * Compiling with -DNSLACK disable this feature, prioritizing memory over speed. */
#ifndef NSLACK
	#define NUM_SLACK_POS 1
	#define idx(n) (n)
#else
	#define NUM_SLACK_POS 0
	#define idx(n) ((n)-1)
#endif

#define weight(n) weights[idx(n)] 
#define low_max_inharm(n) low_max_inharms[idx(n)] 
#define high_max_inharm(n) high_max_inharms[idx(n)] 
#define candidate(n) candidates[idx(n)] 
#define series(n) series[idx(n)] 
#define partial(n) partials[idx(n)] 

/* PDA-related data structure */

typedef struct partial_t {
	double frq;
	double mag;
} Partial;

typedef enum {PARTIAL, REFERENCE} HarmonicType;

typedef struct harmonic_t {
	HarmonicType type;
	union part_ref_t {
		Partial*  partial;
		Partial** ref;
	} harm;
} Harmonic;

typedef struct candidate_t {
	double         prom; /* prominence (Eq. 8) */
	double         wahm; /* weighted average harmonic magnitude (Eq. 13) */
	unsigned short highest_possible_harm;
	Harmonic*      series;
} Candidate;

/* Complex-related data structure */

typedef struct complex_t {
    double re;
    double im;
} Complex;

/* general constants */
typedef enum {FALSE, TRUE} Boolean;
enum {OK, ERROR};
const static unsigned char MAX_ERR_MSG_LEN = 255;

/* trigonometric constants */
const static double TWO_PI = 6.283185307179586476925286766559;
const static double PI = 3.1415926535897932384626433832795;

enum {MIN_BIT_DEPTH = 8, MAX_BIT_DEPTH = 24};
enum {MIN_FFT_SIZE = 128, MAX_FFT_SIZE = 16384};
enum {MIN_SAMPLE_RATE = 128, MAX_SAMPLE_RATE = 192000};
const static unsigned char MAX_DYNAMIC_RANGE = 144;
const static unsigned short HARM_MYSTIC = 16;
const static double MAX_PEAK_FREQ_ADJUSTMENT = 0.5; /* 0.5 is the obvious choice, but 0.52 is also possible */
const static double THIRD_OCTAVE_UP_LOG10 = 0.10034333188799373173791296490816;
const static double QUARTER_TONE_UP   = 1.02930223664349;
const static double QUARTER_TONE_DOWN = 0.971531941153606;
const static double C8_FREQ = 4186.0090448096; /* freq. (Hz) of grand piano highest note */
const static double MIN_EXPECTED_F0_LOWER_BOUND = 20; /* minimum audible frequency */

/* other constants */
#define NULL_PARTIAL_FREQ 0.123456789

/* Unit Testing functions */
extern void   pda_set_partials(const Partial p[], const unsigned short sz);

/* Debugging functions */
extern void     pda_print_spec(void);
extern double   pda_frame_energy(const int frame[]);
extern double   pda_spec_energy(void);
extern void     pda_print_inharms(void);
extern void     pda_print_weights(void);
extern void     pda_print_partials(void);
extern void     pda_print_candidates(void);
extern Partial* pda_get_strongest(void);

#endif /*_MITRE_PDA_H_*/
