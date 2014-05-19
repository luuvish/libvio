
void sei_rbsp(rbsp_t *rbsp) {
	do {
		sei_message(rbsp);
	} while (more_rbsp_data(rbsp));

	rbsp_trailing_bits(rbsp);
}

void sei_message(rbsp_t *rbsp) {
	payloadType = 0;
	while (rbsp->next_bits(8) == 0xff) {
		uint32_t ff_byte = rbsp->f(8);
		assert(ff_byte == 0xff);

		payloadType += 255;
	}
	last_payload_type_byte = rbsp->u(8);
	payloadType += last_payload_type_byte;

	payloadSize = 0;
	while (rbsp->next_bits(8) == 0xff) {
		uint32_t ff_byte = rbsp->f(8);
		assert(ff_byte == 0xff);
		payloadSize += 255;
	}
	last_payload_size_byte = rbsp->u(8);
	payloadSize += last_payload_size_byte;

	sei_payload(payloadType, payloadSize);
}
