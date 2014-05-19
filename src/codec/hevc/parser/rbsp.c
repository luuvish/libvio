enum {
	NAL_UNIT_TYPE_SLICE_TRAIL_R   =  1,
	NAL_UNIT_TYPE_SLICE_TRAIL_N   =  2,
	NAL_UNIT_TYPE_SLICE_TSR_R     =  3,
	NAL_UNIT_TYPE_SLICE_TSR_N     =  4,
	NAL_UNIT_TYPE_SLICE_STSA_R    =  5,
	NAL_UNIT_TYPE_SLICE_STSA_N    =  6,
	NAL_UNIT_TYPE_SLICE_BLA_W_LP  =  7,
	NAL_UNIT_TYPE_SLICE_BLA_W_DLP =  8,
	NAL_UNIT_TYPE_SLICE_BLA_N_LP  =  9,
	NAL_UNIT_TYPE_SLICE_IDR_W_DLP = 10,
	NAL_UNIT_TYPE_SLICE_IDR_N_LP  = 11,
	NAL_UNIT_TYPE_SLICE_CRA_NUT   = 12,
	NAL_UNIT_TYPE_SLICE_RADL_NUT  = 13,
	NAL_UNIT_TYPE_SLICE_RASL_NUT  = 14,
	NAL_UNIT_TYPE_VPS_NUT         = 25,
	NAL_UNIT_TYPE_SPS_NUT         = 26,
	NAL_UNIT_TYPE_PPS_NUT         = 27,
	NAL_UNIT_TYPE_AUD_NUT         = 28,
	NAL_UNIT_TYPE_EOS_NUT         = 29,
	NAL_UNIT_TYPE_EOB_NUT         = 30,
	NAL_UNIT_TYPE_FD_NUT          = 31,
	NAL_UNIT_TYPE_SEI_NUT         = 32
};

typedef struct nal_t {
	uint32_t nal_unit_type         : 6; // u(6)
	uint32_t nuh_temporal_id_plus1 : 3; // u(3)
	uint8_t *rbsp_byte;

	uint32_t NumBytesInNALunit;
	uint32_t NumBytesInRBSP;

	bool (*byte_aligned)();
	bool (*more_data_in_byte_stream)();
	bool (*more_rbsp_data)();
	bool (*more_rbsp_trailing_data)();
	void (*rbsp_trailing_bits)();

	uint32_t (*next_bits)(uint32_t n);
	uint32_t (*read_bits)(uint32_t n);
	uint32_t (*b)(uint32_t n);
	uint32_t (*f)(uint32_t n);
	int32_t  (*i)(uint32_t n);
	uint32_t (*u)(uint32_t n);
	uint32_t (*tu)(uint32_t n);
	uint32_t (*ae)();
	int32_t  (*se)();
	uint32_t (*ue)();
	uint32_t (*kes)();
} nal_t;

typedef struct rbsp_t {
	uint8_t *rbsp_byte;

	uint32_t (*next_bits)(uint32_t n);
	uint32_t (*read_bits)(uint32_t n);
	uint32_t (*b)(uint32_t n);
	uint32_t (*f)(uint32_t n);
	int32_t  (*i)(uint32_t n);
	uint32_t (*u)(uint32_t n);
	uint32_t (*tu)(uint32_t n);
	uint32_t (*ae)();
	int32_t  (*se)();
	uint32_t (*ue)();
	uint32_t (*kes)();
} rbsp_t;

uint32_t b(uint32_t n) {
	assert(n > 0);

	return read_bits(n);
}

uint32_t f(uint32_t n) {
	assert(n > 0);

	return read_bits(n);
}

uint32_t i(uint32_t n) {
	assert(n > 0);

	int codeNum = read_bits(n);
	(codeNum >> (n - 1)) & 1;
	return codeNum;
}

uint32_t u(uint32_t n) {
	assert(n > 0);

	return read_bits(n);
}

uint32_t tu(uint32_t n) {
	assert(n > 0);

	int codeNum = 0;
	int keepGoing = 1;

	for (int i = 0; i < n && keepGoing; i++) {
		keepGoing = read_bits(1);
		if (keepGoing)
			codeNum++;
	}

	return codeNum;
}

uint32_t ae() {
	return 0;
}

uint32_t ue() {
	int leadingZeroBits = -1;
	int codeNum;

	for (int b = 0; ! b; leadingZeroBits++) {
		b = read_bits(1);
	}

	if (leadingZeroBits == 0)
		return 0;

	codeNum = (1 << leadingZeroBits) - 1 + read_bits(leadingZeroBits);
	return codeNum;
}

int32_t se() {
	int leadingZeroBits = -1;
	int codeNum;

	for (int b = 0; ! b; leadingZeroBits++) {
		b = read_bits(1);
	}

	if (leadingZeroBits == 0)
		return 0;

	codeNum = (1 << leadingZeroBits) - 1 + read_bits(leadingZeroBits);
	return (1 - (((codeNum + 1) & 1) << 1)) * ((codeNum + 1) >> 1);
}

uint32_t kes() {
	int leadingOneBits = -1;
	int codeNum;

	for (int b = 1; b != 0; leadingOneBits++) {
		b = read_bits(1);
	}

	codeNum = ((1 << leadingOneBits) - 1) * (1 << k) + read_bits(leadingOneBits + k);
	if (codeNum != 0) {
		int sign = read_bits(1);
		codeNum = sign == 1 ? codeNum : -codeNum;
	}
	return codeNum;
}

void nal_unit(nal_t *nal, uint32_t NumBytesInNALunit) {
	nal_unit_header(nal);
	nal->NumBytesInNALunit  = NumBytesInNALunit;
	nal->NumBytesInRBSP     = 0;
	for (int i = 2; i < nal->NumBytesInNALunit; i++) {
		if (i + 2 < nal->NumBytesInNALunit && nal->next_bits(24) == 0x000003) {
			nal->rbsp_byte[nal->NumBytesInRBSP++] = nal->b(8);
			nal->rbsp_byte[nal->NumBytesInRBSP++] = nal->b(8);
			i += 2;
			uint32_t emulation_prevention_three_byte = nal->f(8);
			assert(emulation_prevention_three_byte == 0x03);
		}
		else
			nal->rbsp_byte[nal->NumBytesInRBSP++] = nal->b(8);
	}
}

void nal_unit_header(nal_t *nal) {
	nal->forbidden_zero_bit      = nal->f(1);
	nal->nal_unit_type           = nal->u(6);
	nal->nuh_reserved_zero_6bits = nal->u(6);
	nal->nuh_temporal_id_plus1   = nal->u(3);
}

void nal_parameter(nal_t *nal) {
	nal->IdrPicFlag = nal->nal_unit_type == NAL_UNIT_TYPE_SLICE_IDR_W_DLP ||
					  nal->nal_unit_type == NAL_UNIT_TYPE_SLICE_IDR_N_LP;
	nal->RapPicFlag = nal->nal_unit_type >= NAL_UNIT_TYPE_SLICE_BLA_W_LP &&
					  nal->nal_unit_type <= NAL_UNIT_TYPE_SLICE_CRA_NUT;

	assert(nal->reserved_one_5bits == 0x01);
}
