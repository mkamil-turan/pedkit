# Contributing to pedkit

Bug reports, questions and pull requests are welcome on the
[issue tracker](https://github.com/mkamil-turan/pedkit/issues).

## Reporting a problem

Please include the pedkit version (`pedkit --version`), the Python version,
the command you ran, and a minimal `.pedx` file that shows the problem. Use
invented example families only: do not attach real patient data.

## Development setup

```sh
git clone https://github.com/mkamil-turan/pedkit.git
cd pedkit
python -m venv .venv
.venv/bin/pip install -e ".[all,test]"
.venv/bin/python -m pytest -q
```

## Pull requests

- Keep the core (`pedkit.io`, `pedkit.validate`, `pedkit.vis.svg`) free of
  third-party dependencies; numpy, matplotlib and PyQt5 are optional extras.
- Keep the SVG output within SVG 1.2 Tiny (no clip paths, patterns or markers).
- A change to a computation needs a test that compares it with an independent
  result (a hand derivation, the other exact method, simulation or numerical
  integration); see `docs/validation.md`.
- Reports print every result; they do not apply decision thresholds.
- Add an entry to `CHANGELOG.md`.

## Code of conduct

Be respectful and constructive. Harassment of any kind is not tolerated.
