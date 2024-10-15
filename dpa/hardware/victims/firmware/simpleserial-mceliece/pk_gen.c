#include <string.h>

#include "pk_gen.h"

unsigned char mat[ PK_NROWS ][ SYS_N/8 ];

// https://gist.github.com/tijnkooijmans/10981093
uint16_t crc16(char* pData, int length)
{
    uint8_t i;
    uint16_t wCrc = 0xffff;
    while (length--) {
        wCrc ^= *(unsigned char *)pData++ << 8;
        for (i=0; i < 8; i++)
            wCrc = wCrc & 0x8000 ? (wCrc << 1) ^ 0x1021 : wCrc << 1;
    }
    return wCrc & 0xffff;
}

uint16_t chkmat()
{
    return crc16((char*) mat, sizeof(mat));
}

uint16_t strmat(uint16_t off, uint16_t len, uint8_t* data)
{
    uint8_t* matbegin = (uint8_t*) mat;
    // In case we have more rows than space, clamp.
    if (off + len > sizeof(mat))
        len = sizeof(mat) - off;
    memcpy(matbegin + off, data, len);
    return len;
}

uint16_t lodmat(uint16_t off, uint16_t len, uint8_t* data)
{
    uint8_t* matbegin = (uint8_t*) mat;
    // In case we have more rows than space, clamp.
    if (off + len > sizeof(mat))
        len = sizeof(mat) - off;
    memcpy(data, matbegin + off, len);
    return len;
}

// Original code follows here.

// ./subroutines/crypto_declassify.h
static void crypto_declassify(void *x,unsigned long long n)
{
}

// ./subroutines/crypto_uint64.h
typedef uint64_t crypto_uint64;

typedef int64_t crypto_uint64_signed;

static crypto_uint64_signed crypto_uint64_signed_negative_mask(crypto_uint64_signed crypto_uint64_signed_x)
{
	return crypto_uint64_signed_x >> 63;
}

static crypto_uint64 crypto_uint64_nonzero_mask(crypto_uint64 crypto_uint64_x)
{
	return crypto_uint64_signed_negative_mask(crypto_uint64_x) | crypto_uint64_signed_negative_mask(-crypto_uint64_x);
}

static crypto_uint64 crypto_uint64_zero_mask(crypto_uint64 crypto_uint64_x)
{
	return ~crypto_uint64_nonzero_mask(crypto_uint64_x);
}

static crypto_uint64 uint64_is_zero_declassify(uint64_t t)
{
	crypto_uint64 mask = crypto_uint64_zero_mask(t);
	crypto_declassify(&mask,sizeof mask);
	return mask;
}

// pk_gen.c

/* input: secret key sk */
/* output: public key pk */
int ge(void)
{
        int i, j, k;
        int row, c;

        unsigned char mask;
        unsigned char b;


        // gaussian elimination

        for (i = 0; i < (PK_NROWS + 7) / 8; i++)
        for (j = 0; j < 8; j++)
        {
                row = i*8 + j;

#if 1
// We only run the first t+1 columns, after that we can use support splitting.
#define MAX_COLUMNS (SYS_T+1)
                if (row >= MAX_COLUMNS)
                        break;
#else
                if (row >= PK_NROWS)
                        break;
#endif

                for (k = row + 1; k < PK_NROWS; k++)
                {
                        mask = mat[ row ][ i ] ^ mat[ k ][ i ];
                        mask >>= j;
                        mask &= 1;
                        mask = -mask;

                        for (c = 0; c < SYS_N/8; c++)
                                mat[ row ][ c ] ^= mat[ k ][ c ] & mask;
                }

                if ( uint64_is_zero_declassify((mat[ row ][ i ] >> j) & 1) ) // return if not systematic
                {
                        return -1;
                }

                for (k = 0; k < PK_NROWS; k++)
                {
                        if (k != row)
                        {
                                mask = mat[ k ][ i ] >> j;
                                mask &= 1;
                                mask = -mask;

                                for (c = 0; c < SYS_N/8; c++)
                                        mat[ k ][ c ] ^= mat[ row ][ c ] & mask;
                        }
                }
        }

	return 0;
}
