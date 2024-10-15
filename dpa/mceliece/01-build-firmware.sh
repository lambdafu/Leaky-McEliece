#!/bin/bash

# Can be invoked with:
# VERBOSE=1 (show build commands, the default uses -Os)
# clean (clean build files)

make -C ../hardware/victims/firmware/simpleserial-mceliece PLATFORM="CW308_STM32L5" CRYPTO_TARGET=NONE "$@"

