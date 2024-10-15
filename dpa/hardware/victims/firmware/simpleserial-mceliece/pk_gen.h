/*
  This file is for public-key generation
*/

#ifndef PK_GEN_H
#define PK_GEN_H

#include <stdint.h>

#define PK_NROWS (128*13)
// We only run the first 129 columns, rounded up to multiple of 8.
#define SYS_N (128+8)
#define SYS_T 128

// Return checksum of current matrix.
uint16_t chkmat();

// Write LEN bytes from data to matrix at offset OFF. Return bytes written.
uint16_t strmat(uint16_t off, uint16_t len, uint8_t* data);

// Read LEN bytes from matrix to data at offset OFF. Return bytes read.
uint16_t lodmat(uint16_t off, uint16_t len, uint8_t* data);

// Run gaussian elimination.
int ge(void);

#endif
