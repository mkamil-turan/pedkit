"""PyQt5 pedigree viewer (QSvgWidget), with a multifactorial risk tab."""

from __future__ import annotations

import argparse
import sys

from ..io import read
from ..model import Level, Pedigree
from ..validate import validate
from .svg import draw


def _load_qt():
    try:
        from PyQt5 import QtCore, QtSvg, QtWidgets
    except ImportError as e:
        raise SystemExit(
            f"Could not load PyQt5/QtSvg ({e}).\n"
            "Install with: pip install 'pedkit[qt]' (on Linux the system Qt SVG "
            "library may also be needed, e.g. qt5-svg)") from e
    return QtCore, QtSvg, QtWidgets


def _mf_panel(ped: Pedigree, QtCore, QtWidgets):
    """Risk tab when the file has a mode=MF disease; none without numpy/matplotlib."""
    from .mf_panel import mf_diseases
    if not mf_diseases(ped):
        return None
    try:
        from .mf_panel import create_mf_panel
        return create_mf_panel(ped, QtCore, QtWidgets)
    except ImportError as e:
        print(f"Could not open the multifactorial risk tab ({e}); numpy and matplotlib "
              "are required.", file=sys.stderr)
        return None


def create_window(ped: Pedigree, title: str = "pedkit"):
    QtCore, QtSvg, QtWidgets = _load_qt()

    class PedigreeWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self.ped = ped
            self.setWindowTitle(title)
            self.svg = QtSvg.QSvgWidget()
            self.svg.load(QtCore.QByteArray(draw(ped).encode("utf-8")))
            size = self.svg.renderer().defaultSize()
            self.svg.setFixedSize(size)
            scroll = QtWidgets.QScrollArea()
            scroll.setWidget(self.svg)
            scroll.setAlignment(QtCore.Qt.AlignCenter)
            scroll.setStyleSheet("background: white;")
            self.mf_panel = _mf_panel(ped, QtCore, QtWidgets)
            if self.mf_panel is None:
                self.setCentralWidget(scroll)
                self.resize(min(size.width() + 40, 1400), min(size.height() + 80, 900))
            else:
                tabs = QtWidgets.QTabWidget()
                tabs.addTab(scroll, "Pedigree")
                tabs.addTab(self.mf_panel, "Multifactorial risk")
                self.setCentralWidget(tabs)
                self.tabs = tabs
                self.resize(1440, 880)

            menu = self.menuBar().addMenu("&File")
            save_action = menu.addAction("Save &as (SVG/PDF/PNG)…")
            save_action.setShortcut("Ctrl+S")
            save_action.triggered.connect(self.save)
            quit_action = menu.addAction("&Quit")
            quit_action.setShortcut("Ctrl+Q")
            quit_action.triggered.connect(self.close)

            findings = validate(ped)
            n_err = sum(f.level is Level.ERROR for f in findings)
            n_warn = len(findings) - n_err
            self.statusBar().showMessage(
                f"{len(ped.individuals)} individuals · {n_err} errors · {n_warn} warnings")
            if findings:
                self.statusBar().setToolTip("\n".join(str(f) for f in findings))

        def save(self):
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save pedigree", "pedigree.svg", "SVG (*.svg);;PDF (*.pdf);;PNG (*.png)")
            if not path:
                return
            from .export import save
            try:
                save(self.ped, path)
            except (OSError, RuntimeError, ValueError) as e:
                QtWidgets.QMessageBox.warning(self, "Could not save", str(e))

    return PedigreeWindow()


def show(path: str) -> int:
    """Read a file, print the findings to the terminal and show it in a window."""
    _, _, QtWidgets = _load_qt()
    ped = read(path)
    for f in validate(ped):
        print(f, file=sys.stderr)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    window = create_window(ped, f"pedkit — {path}")
    window.show()
    return app.exec_()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="pedkit-view", description="Show a PEDX pedigree.")
    ap.add_argument("file", help=".pedx file")
    a = ap.parse_args(argv)
    return show(a.file)


if __name__ == "__main__":
    sys.exit(main())
