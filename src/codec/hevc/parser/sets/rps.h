#ifndef __RPS_H__
#define __RPS_H__

#define MAX_SHORT_TERM_REF_PIC_SETS (64)
#define MAX_LONG_TERM_REF_PIC_SETS  (256)

typedef struct rps_t {
	uint32_t numDeltaPocs;
	uint32_t numNegativePics;
	uint32_t numPositivePics;

	int32_t  deltaPoc       [MAX_SHORT_TERM_REF_PIC_SETS];
	int32_t  deltaPocS0     [MAX_SHORT_TERM_REF_PIC_SETS];
	int32_t  deltaPocS1     [MAX_SHORT_TERM_REF_PIC_SETS];
	uint32_t usedByCurrPicS0[MAX_SHORT_TERM_REF_PIC_SETS];
	uint32_t usedByCurrPicS1[MAX_SHORT_TERM_REF_PIC_SETS];

} rps_t;

#endif /* __RPS_H__ */
