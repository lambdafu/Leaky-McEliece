echo parameter,trial,seconds 1>&2
# We set the error to 0 to guarantee success. The runtime does not depend on the error for this outcome.
for nt in 51220 102450 348864 460896 6688128 6960119 8192128; do
  for m in botan mceliece; do
    for t in {0..9}; do
      echo -n $m$nt,$t, 1>&2
      (cd $m$nt; /usr/bin/time -f '%e'  ../../build/Stockfish --threads=128 recover --error=0 $m$nt-e0.300000-t00$t.leak > /dev/null )
    done
  done
done

