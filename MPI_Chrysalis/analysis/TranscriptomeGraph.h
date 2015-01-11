#ifndef _TRANSCRIPTOMEGRAPH_H_
#define _TRANSCRIPTOMEGRAPH_H_


//#include <string.h>

#include "analysis/DNAVector.h"


int TranscriptomeGraph(vecDNAVector & seq,
		       FILE * pOut,
		       int k,
		       bool connect = true);








#endif 

