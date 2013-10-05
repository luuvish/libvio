#ifndef _REPORT_H_
#define _REPORT_H_


#include <cstdint>
#include <chrono>

struct SNRParameters {

    // Timing related variables
    std::chrono::system_clock::time_point start_time;
    std::chrono::system_clock::time_point end_time;
    int64_t                               tot_time;

    float       snr [3];
    float       snr1[3];
    float       snra[3];
    float       sse [3];
    float       msse[3];

    // B pictures
    int         Bframe_ctr;
    int         frame_ctr;
    int         frame_no;
    int         g_nFrame;

    int         idr_psnr_number;
};


#endif // _REPORT_H_
