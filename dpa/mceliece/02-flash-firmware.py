#!/bin/env python

import chipwhisperer as cw
PLATFORM="CW308_STM32L5"

scope = cw.scope()
scope.default_setup()

cw.program_target(scope, cw.programmers.STM32FProgrammer,
                  "../hardware/victims/firmware/simpleserial-mceliece/simpleserial-mceliece-CW308_STM32L5.hex")

print(scope.errors)
scope.errors.clear()
scope.dis()
