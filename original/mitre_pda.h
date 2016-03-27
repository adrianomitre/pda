#ifndef MITRE_PDA_H_
#define MITRE_PDA_H_

typedef struct pda_parameters_t PDA_Parameters;

struct pda_parameters_t {
	unsigned short frame_len;
	unsigned int   sample_rate;
	unsigned int   bit_depth;
	unsigned short max_submult;     /* maximum submultiple (i.e., maximum n in Eq. 3) */
	double         dynamic_range;   /* magnitude of full amplitude partial, in dB */
	double         min_expected_f0; /* minimum expected F0, in Hertz */
	double         max_expected_f0; /* maximum expected F0, in Hertz */
	double         min_rel_prom;    /* minimum relative prominence (Beta in Eq. 12) */
};

/* default values for PDA_Parameters fields */
const static unsigned short DEFAULT_FRAME_LEN       = 1024;
const static unsigned short DEFAULT_SAMPLE_RATE     = 44100;
const static unsigned short DEFAULT_BIT_DEPTH       = 16;
const static unsigned short DEFAULT_MAX_SUBMULTIPLE = 20;
const static double         DEFAULT_DYNAMIC_RANGE   = 50;
const static double         DEFAULT_MIN_REL_PROM    = 0.97;
/* The default values for MIN/MAX_EXPECTED_F0 are determined in runtime,
 * according to the sample_rate and frame_len */

/* PDA functions */
extern void set_default_values(PDA_Parameters*);
extern unsigned char prepare_pda(const PDA_Parameters*);
extern void setup_dynamic_range(const double /* range */);
extern double pda(const int[] /* frame */);
extern double pda_prom(const int[] /* frame */, double* /* prom */);
extern double pda_wahm(const int[] /* frame */, double* /* wahm */);
extern double pda_prom_wahm(const int[] /* frame */, double* /* prom */, double* /* wahm */);
extern void finalize_pda();

#endif /*MITRE_PDA_H_*/
