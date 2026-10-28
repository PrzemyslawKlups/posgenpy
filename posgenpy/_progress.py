"""A progress bar with one status line under it, for the long-running sweeps.

Plain text only - tqdm's own text bar and IPython's display, no widgets - so a
notebook shows it without fetching any front-end code.
"""
from time import perf_counter

from tqdm import tqdm

BAR_FORMAT = "{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]"


class _Bar(tqdm):
    # no monitor thread: the sweeps fork their worker processes while the bar
    # is open, and a thread holding a lock at fork time can deadlock the child
    monitor_interval = 0


class _DisplayLine:
    """A file for tqdm that shows its latest line in one notebook output.

    tqdm redraws a line by writing "\\r" and the new text; a notebook would
    append every redraw. This keeps one display and replaces its text instead.
    """

    encoding = "utf-8"   # tells tqdm it may draw the bar in block characters

    def __init__(self):
        from IPython.display import Pretty, display
        self._pretty = Pretty
        self._handle = display(Pretty(""), display_id=True)

    def write(self, text):
        line = text.rsplit("\r", 1)[-1].rstrip("\n")
        if line.strip():
            self._handle.update(self._pretty(line))

    def flush(self):
        pass


def _in_notebook():
    """True inside a Jupyter kernel (VSCode, JupyterLab)."""
    try:
        from IPython import get_ipython
    except ImportError:
        return False
    return getattr(get_ipython(), "kernel", None) is not None


class SweepProgress:
    """A text progress bar over `total` steps, and one line for the latest message.

    The line is overwritten rather than appended to, so a sweep that reports
    every cell and every control dataset shows two lines instead of hundreds:
    how far it has got, and what it did last.

    **No overall time estimate.** A sweep's steps cost very different amounts
    - a real-tip cell writes a per-atom table, a control dataset's cells are
    re-thresholds of arrays it already has - so one rate over the whole bar
    predicts nonsense. Instead each phase (`start_phase`) gets its own
    estimate, from the pace of that phase's own finished steps, shown after
    the elapsed time: "real tip 12/63, ~01:05 left".

    Disabled, `log` is a plain `print` (when `verbose`) and the rest does
    nothing - exactly what the sweep did before it had a bar.
    """

    def __init__(self, total, desc, enabled=True, verbose=True):
        self.verbose = verbose
        self.bar = self._line = self._status = None
        self._phase = None
        if not enabled:
            return

        notebook = _in_notebook()
        self.bar = _Bar(
            total=total, desc=desc, bar_format=BAR_FORMAT,
            file=_DisplayLine() if notebook else None,
            ncols=100 if notebook else None,
        )
        if not verbose:
            return

        if notebook:
            self._status = _DisplayLine()
        else:
            self._line = _Bar(
                total=0, position=abs(self.bar.pos) + 1, bar_format="{desc}"
            )

    def log(self, message):
        if not self.verbose:
            return
        if self._status is not None:
            self._status.write(message)
        elif self._line is not None:
            self._line.set_description_str(message)
        else:
            print(message)

    def update(self, n=1):
        """Move the bar on by `n`, outside any phase (work skipped on resume)."""
        if self.bar is not None and n:
            self.bar.update(n)

    def start_phase(self, name, n_steps):
        """Time the next `n_steps` calls to `step` as one phase."""
        self._phase = dict(name=name, total=n_steps, done=0, started=perf_counter())
        self._show_phase()
        if self.bar is not None:
            self.bar.refresh()

    def step(self, units):
        """One step of the current phase done, worth `units` on the bar."""
        if self._phase is not None:
            self._phase["done"] += 1
            self._show_phase()
        self.update(units)

    def _show_phase(self):
        if self.bar is None:
            return
        phase = self._phase
        text = f"{phase['name']} {phase['done']}/{phase['total']}"
        if 0 < phase["done"] < phase["total"]:
            per_step = (perf_counter() - phase["started"]) / phase["done"]
            left = per_step * (phase["total"] - phase["done"])
            text += f", ~{tqdm.format_interval(left)} left"
        self.bar.set_postfix_str(text, refresh=False)

    def close(self):
        if self.bar is not None:
            self.bar.set_postfix_str("", refresh=False)
        # the bar first: closed the other way round, a terminal ends with the
        # last message on the bar's row and the finished bar under it
        for bar in (self.bar, self._line):
            if bar is not None:
                bar.close()
