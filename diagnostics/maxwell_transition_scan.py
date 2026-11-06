"""maxwell_transition_scan.py -- grid diagnostic, not certified: fine scan of r0 where the second hump of Delta
appears (q = 1.001 ... 1.02), in search of interior bridges (it located the case q = 1.001, r0 = 1500 M certified in
../interval_arithmetic/maxwell_angle.py)."""
import sys
from multiprocessing import Pool
from maxwell_bridge_scan import one
if __name__ == "__main__":
    cases = [(q, r0) for q in [1.001, 1.003, 1.01, 1.02]
             for r0 in [200., 300., 400., 500., 700., 1000., 1500., 2000., 3000., 5000.]]
    with Pool(10) as p:
        for line in p.imap(one, cases):
            print(line); sys.stdout.flush()
