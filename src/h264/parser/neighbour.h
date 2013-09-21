#ifndef _NEIGHBOUR_H_
#define _NEIGHBOUR_H_


void CheckAvailabilityOfNeighbors(mb_t* mb);
void CheckAvailabilityOfNeighborsCABAC(mb_t* mb);

void getNeighbour       (mb_t* mb, int xN, int yN, int mb_size[2], PixelPos *pix);
void get4x4Neighbour    (mb_t* mb, int xN, int yN, int mb_size[2], PixelPos *pix);
void get4x4NeighbourBase(mb_t* mb, int xN, int yN, int mb_size[2], PixelPos *pix);
void get_mb_pos         (VideoParameters *p_Vid, int mb_addr, int mb_size[2], short *x, short *y);

void get_neighbors(mb_t* mb, PixelPos *block, int mb_x, int mb_y, int blockshape_x);
void check_dp_neighbors(mb_t* mb);


struct neighbour_t {
	void get_mb2pos(slice_t* slice, int mbAddr, int& xI, int& yI);
	mb_t* get_pos2mb(slice_t* slice, int xP, int yP, int& mbAddr);
	int predict_nnz(mb_t* mb, int pl, int i, int j);
};


extern neighbour_t neighbour;

#endif /* _NEIGHBOUR_H_ */
