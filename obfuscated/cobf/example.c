#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "pda.h"

const static unsigned short FRAME_SIZE = 1024;
const static double TWO_PI = 6.283185307179586476925286766559;
const static double EPSILON = 0.0001;

int main(int argc, char **argv) {
	PDA_Parameters params;
	unsigned int i;
	int frame[FRAME_SIZE];
	double fund_frq;
	
	for (i = 0; i < FRAME_SIZE; ++i) {
		frame[i] = 32767*sin(8.0*i/FRAME_SIZE*TWO_PI);
	}

	set_default_values(&params);
	params.bit_depth = 16;
	params.sample_rate = 44100;
	params.frame_len = FRAME_SIZE;
	
	prepare_pda(&params);
	fund_frq = pda(frame);
	finalize_pda();

	assert(fabs(fund_frq - 8.0*(44100.0/FRAME_SIZE)) <= EPSILON);
	
	return 0;
}
