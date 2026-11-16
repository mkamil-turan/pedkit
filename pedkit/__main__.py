"""``python -m pedkit`` runs the command line interface (see ``pedkit.cli``)."""

import sys

from .cli import main

sys.exit(main())
