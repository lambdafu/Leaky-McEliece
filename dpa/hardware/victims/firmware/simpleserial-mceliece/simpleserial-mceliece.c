#include "hal.h"
#include <stdlib.h>

#include "simpleserial.h"

#include "pk_gen.h"

// maximum number of rows to store
#define STORE_NROWS 4

uint8_t chk_mat(uint8_t* dummy, uint8_t len)
{
    uint16_t crc = chkmat();
    uint8_t crcSerial[2];
    crcSerial[0] = (crc >> 8) & 0xff;
    crcSerial[1] = crc & 0xff;
    simpleserial_put('r', 2, crcSerial);
    return 0x00;
}


uint8_t store_mat(uint8_t* data, uint8_t len)
{
    uint16_t row = (data[0] << 8) + data[1];
    uint16_t dlen = STORE_NROWS*(SYS_N)/8;
    uint8_t result[2];
    uint16_t resultlen;
    resultlen = strmat(row*(SYS_N/8), dlen, &data[2]);
    result[0] = (resultlen >> 8) & 0xff;
    result[1] = resultlen & 0xff;
    simpleserial_put('r', 2, result);
    return 0x00;
}

uint8_t load_mat(uint8_t* data, uint8_t len)
{
    uint16_t row = (data[0] << 8) + data[1];
    uint16_t dlen = STORE_NROWS*(SYS_N)/8;
    uint8_t result[dlen];
    uint16_t resultlen;
    resultlen = lodmat(row*(SYS_N/8), dlen, &result[2]);
    result[0] = (resultlen >> 8) & 0xff;
    result[1] = resultlen & 0xff;
    if (resultlen < dlen)
        memset(result + 2 + resultlen, 0, dlen - resultlen);
    simpleserial_put('r', 2 + resultlen, result);
    return 0x00;
}

uint8_t do_ge(uint8_t* dummy, uint8_t len)
{
    trigger_high();
    (void) ge();
    trigger_low();

    return chk_mat(0, 0);
}

uint8_t reset(uint8_t* x, uint8_t len)
{
    return 0x00;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    simpleserial_init();
    simpleserial_addcmd('g', 0, do_ge);
    simpleserial_addcmd('c', 0, chk_mat);
    simpleserial_addcmd('s', 2+STORE_NROWS*(SYS_N/8), store_mat);
    simpleserial_addcmd('l', 2, load_mat);
    simpleserial_addcmd('x', 0, reset);
    while(1)
        simpleserial_get();
}
