#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "pda.h"
#include "_pda.h"

/* FFT-related variables */

static unsigned short *rev_table;
static double *sin_table, *cos_table, *unpack_sin_table, *unpack_cos_table;
static double *hann;
static double *tapered, *pow_spec;
static unsigned short spec_size;
static Complex *packed, *unpacked;
static double max_peak_power_adjustment;
static double pow_threshold, min_abs_power;

/* PDA-related variables */

static PDA_Parameters params;
static double base_level;
static double *weights, *low_max_inharms, *high_max_inharms;
static double max_exp_f0_bins, min_exp_f0_bins;
static unsigned short harm_size;
static Partial null_partial;

static Partial *partials;
static unsigned short actual_partials;

static Partial *strongest;
static double highest_frq;

static Candidate *candidates;
static unsigned short cand_size, first_active_cand, last_active_cand;

static Boolean properly_prepared = FALSE;

/********************************************************************************************/
/* FATAL_ERROR                                                                              */
/********************************************************************************************/
static unsigned char fatal_error(const char* fun, const char* msg) {
    fprintf(stderr, "%s: %s.\n", fun, msg);
    return ERROR;
}

/********************************************************************************************/
/* FATAL_ERROR_PREP                                                                         */
/********************************************************************************************/
static unsigned char fatal_error_prep(const char* msg) {
    return fatal_error("prepare_pda", msg);
}

/********************************************************************************************/
/* EMALLOC                                                                                  */
/********************************************************************************************/
static void* emalloc(const size_t size) {
    void* p;

    if ( ( p = malloc(size) ) != NULL )
        return p;

    fprintf(stderr, "Error allocating memory. Aborting...\n");
    exit ( ERROR );
}

/********************************************************************************************/
/* INIT_REV_TABLE                                                                           */
/********************************************************************************************/
static void init_rev_table(const unsigned short n) {
    unsigned short k, i, j, half_n = n >> 1, reverse;

    rev_table = (unsigned short*)emalloc(n * sizeof(unsigned short));

    for (k = 0; k < n; k++)
    {
        for (i = 1, j = half_n, reverse = 0; j > 0; i = (i << 1), j = (j >> 1))
            if (k & j)
                reverse += i;
        rev_table[k] = reverse;
    }
}

/********************************************************************************************/
/* FREE_REV_TABLE                                                                           */
/********************************************************************************************/
static void free_rev_table(void) {
    free(rev_table);
}

/********************************************************************************************/
/* INIT_SIN_TABLE                                                                           */
/********************************************************************************************/
static void init_sin_table(const unsigned short n) {
    unsigned short m;

    sin_table = (double*)emalloc((n + 1) * sizeof(double));

    /*
    	for (m = 0; m <= n; m++)
    		sin_table[m] = 0;
    */

    for (m = 2; m <= n; m <<= 1)
        sin_table[m] = sin(TWO_PI / m);

    /*
    	const unsigned int half_n = n/2;
    	for (m = 1; m <= half_n; m++)
    	{
    		double aux = sin_table[m];
    		sin_table[m] = sin_table[m + half_n];
    		sin_table[m + half_n] = aux;
    	}
    */
}

/********************************************************************************************/
/* FREE_SIN_TABLE                                                                           */
/********************************************************************************************/
static void free_sin_table(void) {
    free(sin_table);
}

/********************************************************************************************/
/* INIT_COS_TABLE                                                                           */
/********************************************************************************************/
static void init_cos_table(const unsigned short n) {
    unsigned short m;

    cos_table = (double*)emalloc((n + 1)* sizeof(double));

    /*
    	for (m = 0; m <= n; m++)
    		cos_table[m] = 0;
    */

    for (m = 2; m <= n; m <<= 1)
        cos_table[m] = cos(TWO_PI / m);

    /*
    	const unsigned int half_n = n/2;
    	for (m = 1; m <= half_n; m++)
    	{
    		double aux = sin_table[m];
    		cos_table[m] = cos_table[m + half_n];
    		cos_table[m + half_n] = aux;
    	}
    */
}

/********************************************************************************************/
/* FREE_COS_TABLE                                                                           */
/********************************************************************************************/
static void free_cos_table(void) {
    free(cos_table);
}

/********************************************************************************************/
/* IS_POWER_OF_2                                                                            */
/********************************************************************************************/
static Boolean is_power_of_2(unsigned long n) {
    if (!n) return FALSE;

    while (n > 2) {
        if (n%2) return FALSE;
        n /= 2;
    }
    return TRUE;
}

/********************************************************************************************/
/* INIT_UNPACK_SIN_TABLE                                                                    */
/********************************************************************************************/
static void init_unpack_sin_table(const unsigned short n) {
    const unsigned short half_n = n/2;
    unsigned short k;

    unpack_sin_table = (double*)emalloc((half_n) * sizeof(double));

    for (k = 0; k < half_n; k++)
        unpack_sin_table[k] = sin(TWO_PI * k / n);
}


/********************************************************************************************/
/* FREE_UNPACK_SIN_TABLE                                                                    */
/********************************************************************************************/
static void free_unpack_sin_table(void) {
    free(unpack_sin_table);
}


/********************************************************************************************/
/* INIT_UNPACK_COS_TABLE                                                                    */
/********************************************************************************************/
static void init_unpack_cos_table(const unsigned short n) {
    const unsigned short half_n = n/2;
    unsigned short k;

    unpack_cos_table = (double*)emalloc((half_n) * sizeof(double));

    for (k = 0; k < half_n; k++)
        unpack_cos_table[k] = cos(TWO_PI * k / n);
}


/********************************************************************************************/
/* FREE_UNPACK_COS_TABLE                                                                    */
/********************************************************************************************/
static void free_unpack_cos_table(void) {
    free(unpack_cos_table);
}


/********************************************************************************************/
/* PREPARE_FFT                                                                              */
/********************************************************************************************/
static void prepare_fft(const unsigned short n) {
    init_rev_table(n);
    init_sin_table(n);
    init_cos_table(n);
}

/********************************************************************************************/
/* FINALIZE_FFT                                                                             */
/********************************************************************************************/
static void finalize_fft(void) {
    free_rev_table();
    free_sin_table();
    free_cos_table();
}

/********************************************************************************************/
/* PREPARE_REAL_FFT_POWER                                                                   */
/********************************************************************************************/
static void prepare_real_fft_power() {
    init_unpack_cos_table(params.frame_len);
    init_unpack_sin_table(params.frame_len);

    packed = (Complex*)emalloc(params.frame_len/2 * sizeof(Complex));
    unpacked = (Complex*)emalloc(spec_size * sizeof(Complex));

    prepare_fft(params.frame_len/2);
}


/********************************************************************************************/
/* FINALIZE_REAL_FFT_POWER                                                                  */
/********************************************************************************************/
static void finalize_real_fft_power(void) {
    free_unpack_cos_table();
    free_unpack_sin_table();

    free(packed);
    free(unpacked);

    finalize_fft();
}

/********************************************************************************************/
/* PACK_REAL_IN_COMPLEX                                                                     */
/********************************************************************************************/
/* packs a real array into a half-sized complex array before the DFT */
static void pack_real_in_complex(Complex packed[], const double unpacked[], const unsigned short n) {
    const unsigned short half_n = n/2;
    unsigned short i, j;

    for (i = j = 0; i < half_n; i++, j++) {
        packed[i].re = unpacked[j];
        packed[i].im = unpacked[++j];
    }
}

/********************************************************************************************/
/* UNPACK_REAL_FROM_COMPLEX                                                                 */
/********************************************************************************************/
/* unpacks a complex array into a half-PLUS-ONE-sized  real array after the DFT */
static void unpack_real_from_complex(Complex unpacked[],
                                     const Complex packed[], const unsigned short n) {
    const unsigned short half_n = n/2;
    unsigned short k;

    unpacked[0].re = packed[0].im; /* TODO: verify this! (first bin) */
    unpacked[0].im = 0;
    for (k = 1; k < half_n; k++) {
        const double cosine = unpack_cos_table[k], sine = unpack_sin_table[k];
        const unsigned short half_n_minus_k = half_n - k;

        unpacked[k].re = 0.5 * (
                             packed[k].re * (1 + sine)
                             +
                             packed[half_n_minus_k].re * (1 - sine)
                             +
                             cosine * (packed[k].im + packed[half_n_minus_k].im)
                         );
        unpacked[k].im = 0.5 * (
                             packed[k].im * (1 + sine)
                             +
                             packed[half_n_minus_k].im * (sine - 1)
                             -
                             cosine * (packed[k].re - packed[half_n_minus_k].re)
                         );
    }
    unpacked[half_n].re = packed[0].re - packed[0].im;
    unpacked[half_n].im = 0;
}

/********************************************************************************************/
/* FFT                                                                                      */
/********************************************************************************************/
static void fft(const Complex in[], Complex out[], const unsigned short int n) {
    register Complex w_m, w, t, u;
    register unsigned short int k, m, j;

    for (k = 0; k < n; k++)
    {
        out[rev_table[k]].re = in[k].re;
        out[rev_table[k]].im = in[k].im;
    }

    for (m = 2; m <= n; m = (m << 1))
    {
        w_m.re = cos_table[m];
        w_m.im = sin_table[m];

        for (k = 0; k < n; k += m)
        {
            w.re = 1;
            w.im = 0;

            for (j = 0; j < (m >> 1); j++)
            {
                t.re = (w.re * out[k + j + (m >> 1)].re) - (w.im * out[k + j + (m >> 1)].im);
                t.im = (w.re * out[k + j + (m >> 1)].im) + (w.im * out[k + j + (m >> 1)].re);

                u.re = out[k + j].re;
                u.im = out[k + j].im;

                out[k + j].re = u.re + t.re;
                out[k + j].im = u.im + t.im;

                out[k + j + (m >> 1)].re = u.re - t.re;
                out[k + j + (m >> 1)].im = u.im - t.im;

                t.re = (w.re * w_m.re) - (w.im * w_m.im);
                t.im = (w.re * w_m.im) + (w.im * w_m.re);
                w.re = t.re;
                w.im = t.im;
            }
        }
    }
}

/********************************************************************************************/
/* INIT_HANN_WINDOW                                                                         */
/********************************************************************************************/
static void init_hann() {
    const unsigned short size = params.frame_len;
    unsigned short i;

    for (i = 0; i < size; i++) {
        hann[i] = 0.5 + 0.5*cos(PI*(2*i-size)/size);
    }
}




/********************************************************************************************/
/* HANN_POW_SPEC                                                                            */
/********************************************************************************************/
/* expect f in bins (convertion to radians is made internally) */
static double hann_pow_spec(const double f) {
    double mag;

    if (f == 0) {
        return 0.25;
    } else {
        mag = sin(PI*f) / (TWO_PI * f * (1 -f*f));
        return mag*mag;
    }
}

/********************************************************************************************/
/* GRANDKE                                                                                  */
/********************************************************************************************/
/*
 * TODO: avoid recalculating the sqrt() for the same bin, once as right, then as left
 *       neighbour in the case of consecutive local maxima.
 */
static double grandke(const double pow_left, const double pow_cent, const double pow_right) {
    double a;
    if (pow_left > pow_right) {
        a = sqrt(pow_cent/pow_left);
        return (2*a - 1)/(a + 1) - 1;
    } else {
        a = sqrt(pow_right/pow_cent);
        return (2*a - 1)/(a + 1);
    }
}

/********************************************************************************************/
/* EXTRACT_PARTIALS                                                                         */
/********************************************************************************************/
static void extract_partials() {
    unsigned short i;
    double delta_frq, tmp_mag;

    strongest = NULL;
    actual_partials = 0;
    for (i = 1; i < spec_size-1; i++) {
        if (pow_spec[i] >= pow_threshold     /* If bin is potentially above base level */
            && pow_spec[i] >= pow_spec[i-1]  /* and it is a local maximum              */
            && pow_spec[i] > pow_spec[i+1])
        {
            delta_frq = grandke(pow_spec[i-1], pow_spec[i], pow_spec[i+1]);
            if (delta_frq >= -MAX_PEAK_FREQ_ADJUSTMENT && delta_frq < MAX_PEAK_FREQ_ADJUSTMENT) {
            	tmp_mag = 10*log10(pow_spec[i]/hann_pow_spec(delta_frq)) - base_level;
                assert(hann_pow_spec(delta_frq) <= hann_pow_spec(0));
                assert(hann_pow_spec(delta_frq) >= hann_pow_spec(MAX_PEAK_FREQ_ADJUSTMENT));
                if (tmp_mag > 0) { /* check is indeed above base level */
                    ++actual_partials;
	            	partial(actual_partials).mag = tmp_mag;
                    partial(actual_partials).frq = i + delta_frq;
                    if (!strongest || partial(actual_partials).mag > strongest->mag) {
                        strongest = &partial(actual_partials);
                    }
                }
            }
	i++; /* optimization: two local maxima are never adjacent */        
	}
    }
    if (actual_partials) highest_frq = partial(actual_partials).frq;
}

/********************************************************************************************/
/* APPLY_WINDOW                                                                             */
/********************************************************************************************/
static void apply_window(const int original[]) {
    unsigned short i;
    for (i = 0; i < params.frame_len; i++) {
        tapered[i] = original[i] * hann[i];
    }
}

/********************************************************************************************/
/* COMPUTE_POWER_SPEC                                                                       */
/********************************************************************************************/
static void compute_power_spec(void) {
    unsigned short i;

    pack_real_in_complex(unpacked, tapered, params.frame_len);

    fft(unpacked, packed, params.frame_len/2);

    unpack_real_from_complex(unpacked, packed, params.frame_len);

    for (i = 0; i < spec_size; i++) {
        pow_spec[i] = unpacked[i].re*unpacked[i].re + unpacked[i].im*unpacked[i].im;
    }
    pow_spec[spec_size-1] /= 4;
}

/********************************************************************************************/
/* INIT_HARM_WEIGHTS                                                                        */
/********************************************************************************************/
/* TODO: consider precomputing these values, as well as the FFT twidles, at compile-time */
static void init_harm_weights() {
    unsigned short h;
    double a, b;

    weights = (double*)emalloc((harm_size + NUM_SLACK_POS) * sizeof(double));

    for (h = 1; h <= min(4,harm_size); h++) {
        weight(h) = 1;
    }
    for (h = 5; h <= harm_size; h++) {
        a = sqrt(h/(h-1))*(h-1);
        b = sqrt((h+1)/h)*h;
        weight(h) = log10(b/a)/THIRD_OCTAVE_UP_LOG10;
    }
}

/********************************************************************************************/
/* INIT_MAX_INHARM                                                                          */
/********************************************************************************************/
static void init_max_inharm() {
    register unsigned short h;

    low_max_inharms = (double*)emalloc((harm_size + NUM_SLACK_POS)*sizeof(double));
    high_max_inharms = (double*)emalloc((harm_size + NUM_SLACK_POS)*sizeof(double));

    for (h = 1; h <= min(16,harm_size); ++h) {
        low_max_inharm(h) = h * QUARTER_TONE_DOWN;
        high_max_inharm(h) = h * QUARTER_TONE_UP;
    }
    if (harm_size >= 17) {
        low_max_inharm(17) = 17 * QUARTER_TONE_DOWN;
        high_max_inharm(17) = 17 * sqrt((h+1.0)/h);
    }
    for (h = 18; h <= harm_size; ++h) {
        low_max_inharm(h) = h * sqrt((h-1.0)/h);
        high_max_inharm(h) = h * sqrt((h+1.0)/h);
    }
}

/********************************************************************************************/
/* INIT_CANDIDATES                                                                          */
/********************************************************************************************/
static void init_candidates() {
    unsigned short i, j, k;
    
    cand_size = min(floor(params.sample_rate/(2*params.min_expected_f0)), params.max_submult);
    candidates = (Candidate*)emalloc((cand_size + NUM_SLACK_POS)*sizeof(Candidate));

    for (i = 1; i <= cand_size; ++i) {
        candidate(i).series = (Harmonic*)emalloc((harm_size + NUM_SLACK_POS)*sizeof(Harmonic));
        for (j = 1; j <= harm_size; ++j) {
            candidate(i).series(j).type = PARTIAL;
        }
    }

    for (i = 1; i <= min(cand_size, params.max_submult/2); ++i) {
        for (j = 1; j <= HARM_MYSTIC/2; ++j) {
            for (k = 2; k*i <= cand_size && k*j <= min(harm_size, HARM_MYSTIC); ++k) {
            	if (candidate(i*k).series(j*k).type != REFERENCE) { /* prevent chained references */
	                candidate(i*k).series(j*k).type = REFERENCE;
	                candidate(i*k).series(j*k).harm.ref = &candidate(i).series(j).harm.partial;
            	}
            }
        }
    }
}

/********************************************************************************************/
/* FINALIZE_CANDIDATES                                                                      */
/********************************************************************************************/
static void finalize_candidates() {
	unsigned short i;
	
	for (i = 1; i <= cand_size; ++i) {
		free(candidate(i).series);
	}
	free(candidates);
}

/********************************************************************************************/
/* GET_HARM                                                                                 */
/********************************************************************************************/
static Partial* get_harm(const unsigned short c, const unsigned short h) {
    if (candidate(c).series(h).type == REFERENCE) {
        return *(candidate(c).series(h).harm.ref);
    } else {
        return candidate(c).series(h).harm.partial;
    }
}

/********************************************************************************************/
/* CREATE_ACTUAL_CANDIDATES                                                                 */
/********************************************************************************************/
static void create_actual_candidates() {
    unsigned short i, j;
    double cand_frq;
    
    if (actual_partials == 0) {
    	last_active_cand = 0;
        return;
    }
    assert(strongest != NULL);
    
    /* first_active_cand HAS to be calculated before last_active_cand */
    first_active_cand = ceil(strongest->frq / max_exp_f0_bins);
    last_active_cand = min(params.max_submult, floor(strongest->frq / min_exp_f0_bins));
    if (last_active_cand < first_active_cand) {
		last_active_cand = 0;
    	return;
    }
    
	assert(first_active_cand > 0);
	assert(last_active_cand <= harm_size);
	assert(first_active_cand <= last_active_cand);
    
	/* the following loop must start in 1, i.e., include inactive candidates,
	 * because there may be referenced from active ones. */
    for (i = 1; i <= last_active_cand; ++i) {
    	cand_frq = strongest->frq / i;
        candidate(i).highest_possible_harm = log_round(highest_frq / cand_frq); 
        for (j = 1; j <= candidate(i).highest_possible_harm; ++j) {
            if (candidate(i).series(j).type == PARTIAL) {
                candidate(i).series(j).harm.partial = &null_partial;
            }
        	assert( candidate(i).series(j).type == PARTIAL || get_harm(i,j)->mag == 0 );
        	assert( candidate(i).series(j).type == PARTIAL || get_harm(i,j)->frq == NULL_PARTIAL_FREQ );
        }
    }
}

/********************************************************************************************/
/* REL_INHARM                                                                               */
/********************************************************************************************/
static double rel_inharm(const double measured, const double ref) {
	return max(measured,ref) / min(measured,ref);
}

/********************************************************************************************/
/* UPDATE_CAND_HARM                                                                         */
/********************************************************************************************/
static void update_cand_harm(const unsigned short c, const unsigned short p, const unsigned short h) {
    if (partial(p).mag > candidate(c).series(h).harm.partial->mag
        || (partial(p).mag == candidate(c).series(h).harm.partial->mag
        	&& rel_inharm(partial(p).frq, strongest->frq/c*h)
               < rel_inharm(candidate(c).series(h).harm.partial->frq, strongest->frq/c*h) ) )
    {
        candidate(c).series(h).harm.partial = &partial(p);
    }
}

/********************************************************************************************/
/* IS_INHARM_ACCEPTABLE                                                                     */
/********************************************************************************************/
static Boolean is_inharm_acceptable(const double ratio, const unsigned short h) {
    if (ratio > low_max_inharm(h) && ratio < high_max_inharm(h)) {
        return TRUE;
    } else {
        return FALSE;
    }
}

/********************************************************************************************/
/* COLLECT_SERIES                                                                           */
/********************************************************************************************/
static void collect_series() {
    unsigned short i, j,  h;
    double ratio, cand_grs_period;

    /* the loop must begin from 1, NOT from "first_active_cand" because of references */
    for (i = 1; i <= last_active_cand; ++i) {
        cand_grs_period = (double)i / strongest->frq;
        for (j = 1; j <= actual_partials; ++j) {
            assert(partial(j).mag > 0);
            ratio = partial(j).frq * cand_grs_period;
            h = log_round(ratio);
            assert(h <= harm_size);
            if (candidate(i).series(h).type == PARTIAL) {
                if (is_inharm_acceptable(ratio, h)) {
                    update_cand_harm(i, j, h);
                }
            }
        }
    }
}

/********************************************************************************************/
/* FAST_COLLECT_SERIES                                                                      */
/********************************************************************************************/
/* TODO: it seems to be buggy, fix it and USE IT */
static void fast_collect_series() {
    unsigned short i, j, k, h;
    double ratio, cand_grs_period;

    for (i = first_active_cand; i <= last_active_cand; ++i) {
        cand_grs_period = (double)i / strongest->frq;
        for (j = actual_partials; j > 0 ; --j) {
            assert(partial(j).mag > 0);
            ratio = partial(j).frq * cand_grs_period;
            h = log_round(ratio);
            assert(h <= harm_size);
            if (candidate(i).series(h).type == PARTIAL) { /* TODO: indeed necessary? */
                if (is_inharm_acceptable(ratio, h)) {
                    update_cand_harm(i, j, h);
                    for (k = (HARM_MYSTIC/h)+1; k*i <= last_active_cand; ++k) {
                        if (is_inharm_acceptable(k*ratio, k*h)) {
                            update_cand_harm(k*i, j, k*h);
                        }
                    }
                }
            }
        }
    }
}

/********************************************************************************************/
/* COMPUTE_PROMS_AND_WAHMS                                                                  */
/********************************************************************************************/
static void compute_proms_and_wahms() {
    unsigned short i, j;
    double total, tot_weight, aux;
    
    for (i = first_active_cand; i <= last_active_cand; ++i) {
        total = tot_weight = aux = 0;
        for (j = 1; j <= candidate(i).highest_possible_harm; ++j) {
        	if (get_harm(i, j)->mag > 0) {
	            total += get_harm(i, j)->mag * weight(j);
	            tot_weight += weight(j) + aux;
	            aux = 0;
        	} else {
        		aux += weight(j);
        	}
        }
        candidate(i).prom = total;
        candidate(i).wahm = tot_weight > 0 ? total / tot_weight : 0;
    }
}

/********************************************************************************************/
/* REFINED_CAND_FREQ                                                                        */
/********************************************************************************************/
static double refined_cand_freq(const unsigned short c) {
    double tot_weight = 0, ref_frq = 0, curr_weighted_mag;
    unsigned short h;
    Partial* part;

    assert(c >= first_active_cand);
    assert(c <= last_active_cand);
    assert(actual_partials > 0);
    assert(strongest != NULL);
    assert(c <= last_active_cand);
    assert(last_active_cand > 0);
    for (h = 1; h <= candidate(c).highest_possible_harm; ++h) {
        part = get_harm(c, h);
        if (part == &null_partial) {
        	continue;
        }
        assert(part->mag != 0.0 && part->frq != NULL_PARTIAL_FREQ);
        assert(part != NULL);
        curr_weighted_mag = part->mag * weight(h);
        assert(curr_weighted_mag >= 0);
        ref_frq += (part->frq / h) * curr_weighted_mag;
        assert(is_inharm_acceptable(part->frq / (strongest->frq/c), h));
        tot_weight += curr_weighted_mag;
    }
    assert(tot_weight != 0);
    assert(ref_frq/tot_weight > (strongest->frq/c) * QUARTER_TONE_DOWN);
    assert(ref_frq/tot_weight < (strongest->frq/c) * QUARTER_TONE_UP);
    return ref_frq / tot_weight;
}












































/*==========================================================================================*/
/* DEBUGGING                                                                                */
/*==========================================================================================*/

void pda_print_spec() {
    unsigned short i;

    for (i = 0; i <= spec_size; i++) {
        fprintf(stderr, "pow_spec[%d] = %g\n", i, pow_spec[i]);
    }
}

double pda_frame_energy(const int frame[]) {
    double tot_energy = 0;
    unsigned short i;

    for (i = 0; i < params.frame_len; i++) {
        tot_energy += (frame[i]*frame[i]);
    }
    return tot_energy;
}

double pda_spec_energy() {
    double tot_energy = 0;
    unsigned short i;

    for (i = 0; i <= spec_size; i++) {
        tot_energy += pow_spec[i];
    }
    return tot_energy;
}



static void print_inharm(unsigned n) {
    fprintf(stderr, "inharm(%d)= [%g, %g]\n", n, low_max_inharm(n), high_max_inharm(n));
}

void pda_print_inharms() {
    unsigned short i;
    for (i = 1; i <= harm_size; ++i) {
        print_inharm(i);
    }
}

static void print_weight(unsigned n) {
    fprintf(stderr, "w(%d)= %g\n", n, weights[n-1]);
}

void pda_print_weights() {
    unsigned short i;
    for (i = 1; i <= harm_size; ++i) {
        print_weight(i);
    }
}

Partial* pda_get_strongest() {
	return strongest;
}

static void find_strongest() {
    unsigned short i;

    strongest = NULL;
    for (i = 1; i <= actual_partials; ++i) {
        if (!strongest || partial(i).mag > strongest->mag) {
            strongest = &partial(i);
        }
    }
}

static void print_partial(unsigned short i) {
    if (strongest == &partials[i]) {
        fprintf(stderr, "p(%d)= [f: %g, %g] <=== strongest\n", i+1, partials[i].frq, partials[i].mag);
    } else {
        fprintf(stderr, "p(%d)= [f: %g, %g]\n", i+1, partials[i].frq, partials[i].mag);
    }
}


void pda_print_partials(void) {
    unsigned short i;

    for (i = 0; i < actual_partials; i++) {
        print_partial(i);
    }
}

static Boolean is_harm_ref(const unsigned short c, const unsigned short h) {
    return candidate(c).series(h).type == REFERENCE ? TRUE : FALSE;
}

static void print_candidate(unsigned short n) {
    const Candidate* c = &candidate(n);
    Partial* p;
    unsigned short h;

    fprintf(stderr, "cand(%d) => grs_frq: %g, prom: %g, wahm: %g (high.pos.h.%d)\n", n, (strongest->frq)/n, c->prom, c->wahm, c->highest_possible_harm);
    for (h = 1; h <= candidate(n).highest_possible_harm; ++h) {
        p = get_harm(n, h);
        assert(p != NULL);
        if (p->mag != 0) {
            fprintf(stderr, "\th(%d) = [%g, %g] %s\n", h, p->frq, p->mag, (is_harm_ref(n,h) ? "(ref)" : ""));
        }
    }
}

void pda_print_candidates() {
    unsigned short i;

    for (i = first_active_cand; i <= last_active_cand; ++i) {
        print_candidate(i);
    }
}

void pda_set_partials(const Partial p[], const unsigned short sz) {
	unsigned short i;
	assert(sz > 0);
	assert(sz <= spec_size/2);
	assert(p != NULL);
	
	for (i = 1; i <= sz; ++i) {
		partial(i).frq = p[i-1].frq;
		partial(i).mag = p[i-1].mag;
	}
	actual_partials = sz;
	if (actual_partials) {
		highest_frq = partial(actual_partials).frq;
		find_strongest();
	}
}








































/********************************************************************************************/
/* SET_DEFAULT_VALUES                                                                       */
/********************************************************************************************/
void set_default_values(PDA_Parameters *p) {
    /* partials with magnitude 'dynamic_range' dB below maximum will be discarded */
    p->dynamic_range   = DEFAULT_DYNAMIC_RANGE; /* in decibels, s.t. maximum mag. is mapped to 0 dB */
    p->min_rel_prom    = DEFAULT_MIN_REL_PROM;
    p->max_submult = DEFAULT_MAX_SUBMULTIPLE;
    p->frame_len       = DEFAULT_FRAME_LEN;
    p->sample_rate     = DEFAULT_SAMPLE_RATE;
    p->bit_depth       = DEFAULT_BIT_DEPTH;
    p->min_expected_f0 = 0; /* will be set to freq. resolution by prepare_pda  */
    p->max_expected_f0 = 0; /* will be set to Nyquist freq. by prepare_pda     */
}

/********************************************************************************************/
/* SETUP_DYNAMIC_RANGE                                                                       */
/********************************************************************************************/
void setup_dynamic_range(const double range) {
    if (range <= 0 || range > MAX_DYNAMIC_RANGE) {
        fatal_error("pda_dynamic_range", "dynamic range must satisfy 0 < DR <= 144 dB");
    }
    double max_possible_pow = params.frame_len * params.frame_len * pow(4, params.bit_depth-2);
    base_level = 10*log10(max_possible_pow) - range;
    min_abs_power = max_possible_pow * pow(10, -0.1 * range);
}

/********************************************************************************************/
/* PREPARE_PDA                                                                              */
/********************************************************************************************/
/* receives a frame of length at least params. */
unsigned char prepare_pda(const PDA_Parameters* args) {
    const double nyquist = args->sample_rate / 2.0;
    char err_msg[MAX_ERR_MSG_LEN];
    double freq_res;
    
    if (properly_prepared) {
    	return fatal_error_prep("cannot initialize multiple times consecutively");
    }
    
    assert(args != NULL);

    memcpy(&params, args, sizeof(PDA_Parameters));
    
    /* verify parameters consistency */
    if (params.sample_rate < MIN_SAMPLE_RATE || params.sample_rate > MAX_SAMPLE_RATE) {
        sprintf(err_msg,"sample_rate must be between %d and %d", MIN_SAMPLE_RATE, MAX_SAMPLE_RATE);
        return fatal_error_prep(err_msg);
    }
    if (params.frame_len < MIN_FFT_SIZE || params.frame_len > MAX_FFT_SIZE || !is_power_of_2(params.frame_len)) {
        sprintf(err_msg, "frame_len must be a power of 2 within %d and %d", MIN_FFT_SIZE, MAX_FFT_SIZE);
        return fatal_error_prep(err_msg);
    }

    /* compute deductible parameters set to 0 */
    freq_res = 2 * params.sample_rate / (double)params.frame_len;

    if (!params.min_expected_f0) {
        params.min_expected_f0 = freq_res < MIN_EXPECTED_F0_LOWER_BOUND ? MIN_EXPECTED_F0_LOWER_BOUND : freq_res;
    }
    if (!params.max_expected_f0) {
        params.max_expected_f0 = nyquist > C8_FREQ ? C8_FREQ : nyquist;
    }

    /* verify parameters consistency */
    if (params.min_expected_f0 < MIN_EXPECTED_F0_LOWER_BOUND || params.min_expected_f0 > nyquist) {
        sprintf(err_msg, "min_expected_f0 must lie between %g Hz and Nyquist frequency", MIN_EXPECTED_F0_LOWER_BOUND);
        return fatal_error_prep(err_msg);
    }
    if (params.max_expected_f0 < MIN_EXPECTED_F0_LOWER_BOUND || params.max_expected_f0 > nyquist) {
        sprintf(err_msg, "max_expected_f0 must lie between %g Hz and Nyquist frequency", MIN_EXPECTED_F0_LOWER_BOUND);
        return fatal_error_prep(err_msg);
    }
    if (params.min_expected_f0 >= params.max_expected_f0) {
    	return fatal_error_prep("max_expected_f0 must be greater than min_expected_f0");
    }
    if (params.min_expected_f0 < freq_res) {
    	return fatal_error_prep("min_expected_f0 cannot be smaller than the frequency resolution,\n"
                         "which is sample rate divided by half the frame length");
    }
    if (params.dynamic_range <= 0) {
    	return fatal_error_prep("invalid def_norm_factor");
    }
    if (params.bit_depth % 8 || params.bit_depth < MIN_BIT_DEPTH || params.bit_depth > MAX_BIT_DEPTH) {
        sprintf(err_msg, "bit_depth must be a multiple of weight within %d and %d", MIN_BIT_DEPTH, MAX_BIT_DEPTH);
        return fatal_error_prep(err_msg);
    }
    
    /* From here on, every parameter is supposed valid and no error should ocurr. */

    /* expected F0 range in bins, instead of Hertz */
    min_exp_f0_bins = 2*params.min_expected_f0 / freq_res;
    max_exp_f0_bins = 2*params.max_expected_f0 / freq_res;

    /* allocate memory and initialize structures */
    hann = (double*)emalloc(params.frame_len * sizeof(double));
    init_hann();
    tapered = (double*)emalloc(params.frame_len * sizeof(double));
    spec_size = params.frame_len/2 + 1;
    pow_spec = (double*)emalloc(spec_size * sizeof(double));
    partials = (Partial*)emalloc((spec_size/2 + NUM_SLACK_POS)*sizeof(Partial));

    setup_dynamic_range(params.dynamic_range);

    prepare_real_fft_power();

    max_peak_power_adjustment = hann_pow_spec(MAX_PEAK_FREQ_ADJUSTMENT);

    harm_size = log_round(params.sample_rate / (2*params.min_expected_f0));

    init_harm_weights();
    init_max_inharm();

    init_candidates();
    
    null_partial.frq = NULL_PARTIAL_FREQ;
    null_partial.mag = 0.0;

    properly_prepared = TRUE;
    
    return OK;
}

/********************************************************************************************/
/* FINALIZE_PDA                                                                             */
/********************************************************************************************/
void finalize_pda() {
	if (properly_prepared) {
	    finalize_real_fft_power();
	
	    free(partials);
	    free(pow_spec);
	    free(tapered);
	    free(hann);
	
	    free(high_max_inharms);
	    free(low_max_inharms);
	    free(weights);
	    
	    finalize_candidates();
	    
	    properly_prepared = FALSE;
	}
}

/********************************************************************************************/
/* PDA                                                                                      */
/********************************************************************************************/
double pda(const int frame[]) {
    return pda_prom(frame, NULL);
}

/********************************************************************************************/
/* PDA_PROM                                                                                 */
/********************************************************************************************/
double pda_prom(const int frame[], double* prom) {
	return pda_prom_wahm(frame, prom, NULL);
}

/********************************************************************************************/
/* PDA_WAHM                                                                                 */
/********************************************************************************************/
double pda_wahm(const int frame[], double* wahm) {
	return pda_prom_wahm(frame, NULL, wahm);
}

/********************************************************************************************/
/* PDA_PROM_WAHM                                                                            */
/********************************************************************************************/
double pda_prom_wahm(const int frame[], double* prom, double* wahm) {
    unsigned short i, max_prom_cand, winner_cand;
    double min_abs_prom;
    
#ifndef UNIT_TESTING
    assert(frame != NULL);

    apply_window(frame);
    compute_power_spec();
    extract_partials();
#endif /* UNIT_TESTING */
    create_actual_candidates();
    
    if (last_active_cand == 0) { /* no valid candidates */
        if (prom) *prom = 0;
        if (wahm) *wahm = 0;
        return 0;
    }

    collect_series();
    compute_proms_and_wahms();
    
    max_prom_cand = first_active_cand;
    for (i = first_active_cand+1; i <= last_active_cand; ++i) {
        if (candidate(i).prom > candidate(max_prom_cand).prom) {
            max_prom_cand = i;
        }
    }
    min_abs_prom = params.min_rel_prom * candidate(max_prom_cand).prom;
    winner_cand = max_prom_cand;
    for (i = first_active_cand; i <= last_active_cand; ++i) {
        if (candidate(i).prom >= min_abs_prom && candidate(i).wahm > candidate(winner_cand).wahm) {
            winner_cand = i;
        }
    }
    if (prom) *prom = candidate(winner_cand).prom;
    if (wahm) *wahm = candidate(winner_cand).wahm;
    return refined_cand_freq(winner_cand) * params.sample_rate / params.frame_len;
}
