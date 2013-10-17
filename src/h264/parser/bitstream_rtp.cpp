#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

#include <netinet/in.h>

#include "memalloc.h"
#include "bitstream.h"
#include "sets.h"


namespace vio  {
namespace h264 {


#define MAXRTPPAYLOADLEN  (65536 - 40)    //!< Maximum payload size of an RTP packet */
#define MAXRTPPACKETSIZE  (65536 - 28)    //!< Maximum size of an RTP packet incl. header */
#define H264PAYLOADTYPE 105               //!< RTP paylaod type fixed here for simplicity*/
#define H264SSRC 0x12345678               //!< SSRC, chosen to simplify debugging */
#define RTP_TR_TIMESTAMP_MULT 1000        //!< should be something like 27 Mhz / 29.97 Hz */

struct RTPpacket_t {
    uint32_t v;          //!< Version, 2 bits, MUST be 0x2
    uint32_t p;          //!< Padding bit, Padding MUST NOT be used
    uint32_t x;          //!< Extension, MUST be zero
    uint32_t cc;         /*!< CSRC count, normally 0 in the absence
                                  of RTP mixers */
    uint32_t m;          //!< Marker bit
    uint32_t pt;         //!< 7 bits, Payload Type, dynamically established
    uint16_t seq;        /*!< RTP sequence number, incremented by one for
                                  each sent packet */
    uint32_t timestamp;  //!< timestamp, 27 MHz for H.264
    uint32_t ssrc;       //!< Synchronization Source, chosen randomly
    uint8_t* payload;    //!< the payload including payload headers
    uint32_t paylen;     //!< length of payload in bytes
    uint8_t* packet;     //!< complete packet including header and payload
    uint32_t packlen;    //!< length of packet, typically paylen+12
};

int  DecomposeRTPpacket(RTPpacket_t* p);
void DumpRTPHeader(RTPpacket_t* p);
int  RTPReadPacket(RTPpacket_t* p, int bitstream);


void open_rtp(const char *fn, int *p_BitStreamFile)
{
    if (((*p_BitStreamFile) = open(fn, O_RDONLY)) == -1)
        error(500, "Cannot open RTP file '%s'", fn);
}

void close_rtp(int *p_BitStreamFile)
{
    if ((*p_BitStreamFile) != -1) {
        close(*p_BitStreamFile);
        (*p_BitStreamFile) = -1;
    }
}

int get_nalu_from_rtp(nal_unit_t& nal, int BitStreamFile)
{
    static uint16_t first_call = 1;  //!< triggers sequence number initialization on first call
    static uint16_t old_seq = 0;     //!< store the last RTP sequence number for loss detection

    RTPpacket_t* p = new RTPpacket_t;
    p->packet  = new uint8_t[MAXRTPPACKETSIZE];
    p->payload = new uint8_t[MAXRTPPACKETSIZE];

    int ret = RTPReadPacket(p, BitStreamFile);
    nal.forbidden_zero_bit = 1;
    nal.num_bytes_in_nal_unit = 0;

    if (ret > 0) { // we got a packet ( -1=error, 0=end of file )
        if (first_call) {
            first_call = 0;
            old_seq = (uint16_t)(p->seq - 1);
        }

        nal.lost_packets = (uint16_t)(p->seq - (old_seq + 1));
        old_seq = p->seq;

        assert(p->paylen < nal.max_size);

        nal.num_bytes_in_nal_unit = p->paylen;
        memcpy(nal.rbsp_byte, p->payload, p->paylen);
        nal.forbidden_zero_bit = (nal.rbsp_byte[0] >> 7) & 1;
        nal.nal_ref_idc        = (nal.rbsp_byte[0] >> 5) & 3;
        nal.nal_unit_type      = (nal.rbsp_byte[0] & 0x1f);
        if (nal.lost_packets)
            printf("Warning: RTP sequence number discontinuity detected\n");
    }

    // free memory
    delete []p->payload;
    delete []p->packet;
    delete p;

    if (ret < 0)
        nal.num_bytes_in_nal_unit = -1;
    return nal.num_bytes_in_nal_unit;
}

int DecomposeRTPpacket(RTPpacket_t* p)
{
    // consistency check
    assert(p->packlen < 65536 - 28);  // IP, UDP headers
    assert(p->packlen >= 12);         // at least a complete RTP header
    assert(p->payload != NULL);
    assert(p->packet != NULL);

    // Extract header information

    p->v  = (p->packet[0] >> 6) & 0x03;
    p->p  = (p->packet[0] >> 5) & 0x01;
    p->x  = (p->packet[0] >> 4) & 0x01;
    p->cc = (p->packet[0] >> 0) & 0x0F;

    p->m  = (p->packet[1] >> 7) & 0x01;
    p->pt = (p->packet[1] >> 0) & 0x7F;

    memcpy(&p->seq, &p->packet[2], 2);
    p->seq = ntohs((uint16_t)p->seq);

    memcpy(&p->timestamp, &p->packet[4], 4); // change to shifts for unified byte sex
    p->timestamp = ntohl(p->timestamp);
    memcpy(&p->ssrc, &p->packet[8], 4); // change to shifts for unified byte sex
    p->ssrc = ntohl(p->ssrc);

    // header consistency checks
    if (p->v != 2 || p->p != 0 || p->x != 0 || p->cc != 0) {
        printf("DecomposeRTPpacket, RTP header consistency problem, header follows\n");
        DumpRTPHeader(p);
        return -1;
    }
    p->paylen = p->packlen - 12;
    memcpy(p->payload, &p->packet[12], p->paylen);
    return 0;
}

void DumpRTPHeader(RTPpacket_t* p)
{
    for (int i = 0; i < 30; i++)
        printf("%02x ", p->packet[i]);
    printf("Version (V): %d\n", (int)p->v);
    printf("Padding (P): %d\n", (int)p->p);
    printf("Extension (X): %d\n", (int)p->x);
    printf("CSRC count (CC): %d\n", (int)p->cc);
    printf("Marker bit (M): %d\n", (int)p->m);
    printf("Payload Type (PT): %d\n", (int)p->pt);
    printf("Sequence Number: %d\n", (int)p->seq);
    printf("Timestamp: %d\n", (int)p->timestamp);
    printf("SSRC: %d\n", (int)p->ssrc);
}

int RTPReadPacket(RTPpacket_t* p, int bitstream)
{
    int64_t Filepos;
    int     intime;

    assert(p != NULL);
    assert(p->packet != NULL);
    assert(p->payload != NULL);

    Filepos = lseek(bitstream, 0, SEEK_CUR);
    if (4 != read(bitstream, &p->packlen, 4))
        return 0;
    if (4 != read(bitstream, &intime, 4)) {
        lseek(bitstream, Filepos, SEEK_SET);
        printf("RTPReadPacket: File corruption, could not read Timestamp, exit\n");
        exit(-1);
    }

    assert(p->packlen < MAXRTPPACKETSIZE);

    if (p->packlen != (unsigned int)read(bitstream, p->packet, p->packlen)) {
        printf("RTPReadPacket: File corruption, could not read %d bytes\n", (int)p->packlen);
        exit(-1);    // EOF inidication
    }

    if (DecomposeRTPpacket(p) < 0) {
        // this should never happen, hence exit() is ok.  We probably do not want to attempt
        // to decode a packet that obviously wasn't generated by RTP
        printf("Errors reported by DecomposePacket(), exit\n");
        exit(-700);
    }

    assert(p->pt == H264PAYLOADTYPE);
    assert(p->ssrc == H264SSRC);

    return p->packlen;
}


}
}
