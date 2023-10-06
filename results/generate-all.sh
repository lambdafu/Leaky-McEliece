#!/bin/env bash

for kem in botan102450 {mceliece,botan}{348864,460896,6688128,6960119,8192128}; do
  echo "(mkdir -p $kem && cd $kem && sage ../generate-artefacts.sage $kem)"
done

