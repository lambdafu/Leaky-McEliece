#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rng.h"
#include "crypto_kem.h"

unsigned char *leak_buffer = NULL;
size_t leak_buffer_len = 0;
size_t leak_buffer_offset = 0;

void leak_byte(unsigned char out) {
    if (leak_buffer_offset >= leak_buffer_len) {
        leak_buffer_len = leak_buffer_len ? 2*leak_buffer_len : 4096;
        leak_buffer = realloc(leak_buffer, leak_buffer_len);
        assert(leak_buffer);
    }
    leak_buffer[leak_buffer_offset++] = out;
}

void leak_bit(unsigned char bit) {
    leak_byte(bit ? 0x01 : 0x00);
}

void print_bytes(char *prefix, unsigned char *bytes, size_t len) {
    size_t i;
    printf("%s = '", prefix);
    for (i = 0; i < len; i++)
        printf("%02X", bytes[i]);
    printf("'\n");
}

unsigned char hexdigit(char x) {
    if (x >= '0' && x <= '9')
        return x - '0';
    else if (x >= 'a' && x <= 'f')
        return 10 + x - 'a';
    else {
        assert(x >= 'A' && x <= 'F');
        return 10 + x - 'A';
    }
}

unsigned char hex2byte(char *x) {
    return hexdigit(x[0]) * 16 + hexdigit(x[1]);
}

int main(int argc, char* argv[]) {
    int i;
    int ret_val;
    unsigned char seed[48];
    unsigned char *ct = 0;
    unsigned char *ss = 0;
    unsigned char *ss1 = 0;
    unsigned char *pk = 0;
    unsigned char *sk = 0;

    // example seed: '000102030405060708090A0B0C0D0E0F101112131415161718191A1B1C1D1E1F202122232425262728292A2B2C2D2E2F'
    if (argc >= 2) {
        assert(strlen(argv[1]) == 2*48);
        for (i = 0; i < 48; i++)
        seed[i] = hex2byte(&argv[1][2*i]);
    } else {
        for (i = 0; i < 48; i++)
        seed[i] = i;
    }
    randombytes_init(seed, NULL, 256);

    ct = malloc(crypto_kem_CIPHERTEXTBYTES);
    assert(ct);
    ss = malloc(crypto_kem_BYTES);
    assert(ss);
    ss1 = malloc(crypto_kem_BYTES);
    assert(ss1);
    pk = malloc(crypto_kem_PUBLICKEYBYTES);
    assert(pk);
    sk = malloc(crypto_kem_SECRETKEYBYTES);
    assert(sk);
    
    printf("kem = '%s'\n", crypto_kem_PRIMITIVE);
    print_bytes("seed", seed, 48);
    if ( (ret_val = crypto_kem_keypair(pk, sk)) != 0) {
        fprintf(stderr, "crypto_kem_keypair returned <%d>\n", ret_val);
        return -1;
    }
    print_bytes("pk", pk, crypto_kem_PUBLICKEYBYTES);
    print_bytes("sk", sk, crypto_kem_SECRETKEYBYTES);

    if ( (ret_val = crypto_kem_enc(ct, ss, pk)) != 0) {
        fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
        return -1;
    }
    print_bytes("ct", ct, crypto_kem_CIPHERTEXTBYTES);
    print_bytes("ss", ss, crypto_kem_BYTES);

    if ( (ret_val = crypto_kem_dec(ss1, ct, sk)) != 0) {
        fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
        return -1;
    }
    if ( memcmp(ss, ss1, crypto_kem_BYTES) ) {
        fprintf(stderr, "crypto_kem_dec returned bad 'ss' value\n");
        return -1;
    }

    print_bytes("leak", leak_buffer, leak_buffer_offset);
    return 0;
}
