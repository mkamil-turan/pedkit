"""Export to formats other than SVG (PDF, PNG).

Uses the QtSvg module of PyQt5 to avoid another dependency; PyQt5 is
imported only when these functions are called.
"""

from __future__ import annotations

import os

from ..model import Pedigree
from .svg import draw


def _qt():
    try:
        from PyQt5 import QtCore, QtGui, QtSvg
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise RuntimeError(
            "PDF/PNG output needs PyQt5 with QtSvg: pip install 'pedkit[qt]' "
            "(on Linux the system Qt SVG library may also be needed, e.g. qt5-svg)") from e
    if QtGui.QGuiApplication.instance() is None:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        _qt.app = QtGui.QGuiApplication([])  # keep alive for the process lifetime
    return QtCore, QtGui, QtSvg


def _renderer(ped: Pedigree):
    QtCore, QtGui, QtSvg = _qt()
    renderer = QtSvg.QSvgRenderer(QtCore.QByteArray(draw(ped).encode("utf-8")))
    if not renderer.isValid():
        raise RuntimeError("QtSvg could not render the generated SVG.")
    return renderer


def save_pdf(ped: Pedigree, path: str | os.PathLike) -> None:
    QtCore, QtGui, _ = _qt()
    renderer = _renderer(ped)
    size = renderer.defaultSize()
    writer = QtGui.QPdfWriter(os.fspath(path))
    writer.setResolution(72)  # 1 SVG pixel = 1 pt
    writer.setPageSizeMM(QtCore.QSizeF(size.width() * 25.4 / 72, size.height() * 25.4 / 72))
    writer.setPageMargins(QtCore.QMarginsF(0, 0, 0, 0))
    painter = QtGui.QPainter(writer)
    renderer.render(painter, QtCore.QRectF(0, 0, size.width(), size.height()))
    painter.end()


def save_png(ped: Pedigree, path: str | os.PathLike, scale: float = 2.0) -> None:
    QtCore, QtGui, _ = _qt()
    renderer = _renderer(ped)
    size = renderer.defaultSize()
    image = QtGui.QImage(int(size.width() * scale), int(size.height() * scale),
                         QtGui.QImage.Format_ARGB32)
    image.fill(QtCore.Qt.white)
    painter = QtGui.QPainter(image)
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    renderer.render(painter)
    painter.end()
    if not image.save(os.fspath(path)):
        raise RuntimeError(f"Could not save PNG: {path}")


def save(ped: Pedigree, path: str | os.PathLike) -> None:
    """Save as SVG, PDF or PNG depending on the file extension."""
    ext = os.path.splitext(os.fspath(path))[1].lower()
    if ext == ".svg":
        from .svg import save_svg
        save_svg(ped, path)
    elif ext == ".pdf":
        save_pdf(ped, path)
    elif ext == ".png":
        save_png(ped, path)
    else:
        raise ValueError(f"Unsupported extension: {ext!r} (.svg, .pdf, .png)")
