#ifndef _karlin_
#define _karlin_
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
void BlastKarlinBlkCalc(double* scoreProbabilities, int min, int max, float *lambda, float *K, float *H);

int BlastComputeLengthAdjustment(float K,
                             float logK,
                             float alpha_d_lambda,
                             float beta,
                             int query_length,
                             uint32_t db_length,
                             int db_num_seqs,
                             int *length_adjustment);
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif

