"""progress.py -- progress bar for the slow diagnostics, kept out of their logs.

The bar never goes to standard output, so the logs in ../results do not depend on it. If the environment variable
PROGRESS_FILE is set (run_diagnostics.sh sets it), the current line is written to that file, which
`./run_diagnostics.sh --status` prints; otherwise, if standard error is a terminal, the bar is drawn there.
Otherwise nothing is shown. Updates are written at most every two seconds, so the cost is negligible.

    pg = Progress(total, 'label')
    for ...:
        ...
        pg.step()
    pg.close()
"""
import os, sys, time


def _hms(s):
    s = int(round(s)); h, s = divmod(s, 3600); m, s = divmod(s, 60)
    return f'{h}h{m:02d}m' if h else f'{m}m{s:02d}s'


class Progress:
    def __init__(self, total, label=''):
        self.total, self.label, self.done = max(int(total), 1), label, 0
        self.t0 = self.last = time.time()
        self.path = os.environ.get('PROGRESS_FILE')
        self.tty = not self.path and sys.stderr.isatty()
        self._show(force=True)

    def step(self, n=1):
        self.done = min(self.done + n, self.total)
        self._show(force=self.done == self.total)

    def _line(self):
        frac = self.done/self.total; el = time.time() - self.t0
        bar = '#'*int(30*frac) + '-'*(30 - int(30*frac))
        left = f', about {_hms(el/self.done*(self.total - self.done))} left' if self.done else ''
        return f'{self.label} [{bar}] {self.done}/{self.total} ({100*frac:.0f}%), {_hms(el)} elapsed{left}'

    def _show(self, force=False):
        now = time.time()
        if not force and now - self.last < 2: return
        self.last = now
        if self.path:
            tmp = self.path + '.tmp'
            with open(tmp, 'w') as fh: fh.write(self._line() + '\n')
            os.replace(tmp, self.path)
        elif self.tty:
            sys.stderr.write('\r' + self._line()); sys.stderr.flush()

    def close(self):
        self.done = self.total; self._show(force=True)
        if self.tty: sys.stderr.write('\n'); sys.stderr.flush()
