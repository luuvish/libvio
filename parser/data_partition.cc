
#include "global.h"
#include "parser.h"

/*!
 ************************************************************************
 * \brief
 *    Allocates a stand-alone partition structure.  Structure should
 *    be freed by FreePartition();
 *    data structures
 *
 * \par Input:
 *    n: number of partitions in the array
 * \par return
 *    pointer to DataPartition Structure, zero-initialized
 ************************************************************************
 */

DataPartition *AllocPartition(int n)
{
    DataPartition *partArr, *dataPart;
    int i;

    partArr = (DataPartition *)calloc(n, sizeof(DataPartition));
    if (partArr == NULL) {
        snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Data Partition failed");
        error(errortext, 100);
    }

    for (i = 0; i < n; ++i) { // loop over all data partitions
        dataPart = &(partArr[i]);
        dataPart->bitstream = (Bitstream *)calloc(1, sizeof(Bitstream));
        if (dataPart->bitstream == NULL) {
            snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Bitstream failed");
            error(errortext, 100);
        }
        dataPart->bitstream->streamBuffer = (byte *)calloc(MAX_CODED_FRAME_SIZE, sizeof(byte));
        if (dataPart->bitstream->streamBuffer == NULL) {
            snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for streamBuffer failed");
            error(errortext, 100);
        }
    }
    return partArr;
}

/*!
 ************************************************************************
 * \brief
 *    Frees a partition structure (array).
 *
 * \par Input:
 *    Partition to be freed, size of partition Array (Number of Partitions)
 *
 * \par return
 *    None
 *
 * \note
 *    n must be the same as for the corresponding call of AllocPartition
 ************************************************************************
 */
void FreePartition(DataPartition *dp, int n)
{
    int i;

    assert(dp != NULL);
    assert(dp->bitstream != NULL);
    assert(dp->bitstream->streamBuffer != NULL);
    for (i = 0; i < n; ++i) {
        free(dp[i].bitstream->streamBuffer);
        free(dp[i].bitstream);
    }
    free(dp);
}
