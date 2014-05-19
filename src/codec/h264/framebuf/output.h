#ifndef _OUTPUT_H_
#define _OUTPUT_H_


extern void write_stored_frame(VideoParameters *p_Vid, pic_t *fs, int p_out);
extern void direct_output     (VideoParameters *p_Vid, storable_picture *p, int p_out);


#endif //_OUTPUT_H_
