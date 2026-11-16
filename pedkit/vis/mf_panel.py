"""Multifactorial risk tab of the viewer: report and embedded plots."""

from __future__ import annotations

import os

from ..model import Mode, Pedigree


def mf_diseases(ped: Pedigree):
    return [d for d in ped.diseases if d.mode == Mode.MF]


def _default_individual(ped: Pedigree, disease, labels: list[str], phenotypes: list[int]) -> str:
    """First the risk question of the file, then the first individual of unknown phenotype."""
    for q in ped.questions:
        if q.type == "risk" and q.individual and (q.disease in (None, disease.id)):
            for lab in labels:
                if lab == q.individual or lab.endswith("/" + q.individual):
                    return lab
    for lab, f in zip(labels, phenotypes):
        if f == 0:
            return lab
    return labels[0]


def create_mf_panel(ped: Pedigree, QtCore, QtWidgets):
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

    from ..mf.model import prepare_mf
    from ..mf.plots import PLOT_TITLES, figures
    from ..mf.risk import mf_risk

    class Worker(QtCore.QThread):
        """The computation takes a few seconds; it runs in a separate thread
        so the interface does not freeze."""
        done = QtCore.pyqtSignal(object, object)
        failed = QtCore.pyqtSignal(str)

        def __init__(self, options: dict):
            super().__init__()
            self.options = options

        def run(self):
            try:
                r = mf_risk(ped, **self.options)
                self.done.emit(r, figures(r))
            except (KeyError, ValueError, RuntimeError) as e:
                self.failed.emit(str(e))

    class MFPanel(QtWidgets.QWidget):
        def __init__(self):
            super().__init__()
            self.diseases = mf_diseases(ped)
            self.figures: dict = {}
            self.worker = None

            self.disease_box = QtWidgets.QComboBox()
            for d in self.diseases:
                self.disease_box.addItem(f"{d.id}" + (f" · {d.name}" if d.name else ""), d.id)
            self.individual_box = QtWidgets.QComboBox()
            self.method_box = QtWidgets.QComboBox()
            for text, value in (("Elston–Stewart + MCMC", "both"), ("Elston–Stewart", "es"),
                                ("MCMC", "mcmc")):
                self.method_box.addItem(text, value)
            self.h2 = QtWidgets.QDoubleSpinBox()
            self.h2.setRange(0.0, 0.99)
            self.h2.setSingleStep(0.05)
            self.h2.setDecimals(2)
            self.h2.setPrefix("h² = ")
            self.compute_button = QtWidgets.QPushButton("Compute")
            self.save_button = QtWidgets.QPushButton("Save plots…")
            self.save_button.setEnabled(False)
            self.status = QtWidgets.QLabel()

            top = QtWidgets.QHBoxLayout()
            for label, w in (("Disease", self.disease_box), ("Individual", self.individual_box),
                             ("Method", self.method_box)):
                top.addWidget(QtWidgets.QLabel(label))
                top.addWidget(w)
            top.addWidget(self.h2)
            top.addWidget(self.compute_button)
            top.addWidget(self.save_button)
            top.addStretch(1)

            self.report = QtWidgets.QPlainTextEdit()
            self.report.setReadOnly(True)
            self.report.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
            font = QtCore.QCoreApplication.instance().font()
            font.setFamily("monospace")
            font.setStyleHint(font.Monospace)
            font.setPointSizeF(font.pointSizeF() * 0.9)
            self.report.setFont(font)
            self.plot_tabs = QtWidgets.QTabWidget()

            splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
            splitter.addWidget(self.report)
            splitter.addWidget(self.plot_tabs)
            splitter.setStretchFactor(0, 0)       # the plots grow with the window
            splitter.setStretchFactor(1, 1)
            self.splitter = splitter
            # Report lines are at most 76 characters; start with a width that fits them
            self._fit_report_width(["M" * 76])

            layout = QtWidgets.QVBoxLayout(self)
            layout.addLayout(top)
            layout.addWidget(splitter, 1)
            layout.addWidget(self.status)

            self.disease_box.currentIndexChanged.connect(self._disease_changed)
            self.compute_button.clicked.connect(self.compute)
            self.save_button.clicked.connect(self.save)
            self._disease_changed()

        def _fit_report_width(self, lines: list[str]):
            """Widen the left panel so the longest report line fits without scrolling."""
            fm = self.report.fontMetrics()
            text = max(fm.horizontalAdvance(s) for s in lines) if lines else 0
            margin = int(2 * self.report.document().documentMargin()) + 2 * self.report.frameWidth()
            scroll = self.report.verticalScrollBar().sizeHint().width()
            needed = text + margin + scroll + 6
            self.report.setMinimumWidth(needed)
            total = sum(self.splitter.sizes()) or self.width() or 1360
            self.splitter.setSizes([needed, max(total - needed, 1)])

        def _disease(self):
            return self.diseases[self.disease_box.currentIndex()]

        def _disease_changed(self):
            d = self._disease()
            self.individual_box.clear()
            try:
                h = prepare_mf(ped, d)
            except (KeyError, ValueError) as e:
                self.status.setText(f"Could not prepare: {e}")
                self.compute_button.setEnabled(False)
                return
            self.compute_button.setEnabled(True)
            self.individual_box.addItems(h.labels)
            self.individual_box.setCurrentText(_default_individual(ped, d, h.labels, h.phenotype))
            self.h2.setValue(d.heritability if d.heritability is not None else 0.5)
            self.status.setText("" if d.heritability is not None else
                                "No heritability in the file; the h² box value will be used.")

        def compute(self):
            if self.worker is not None and self.worker.isRunning():
                return
            options = dict(individual=self.individual_box.currentText(), disease=self._disease(),
                           method=self.method_box.currentData(),
                           heritability=round(self.h2.value(), 4))
            self.compute_button.setEnabled(False)
            self.status.setText(f"Computing for {options['individual']}…")
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.BusyCursor)
            self.worker = Worker(options)
            self.worker.done.connect(self._done)
            self.worker.failed.connect(self._failed)
            self.worker.finished.connect(self._worker_finished)
            self.worker.start()

        def _worker_finished(self):
            QtWidgets.QApplication.restoreOverrideCursor()
            self.compute_button.setEnabled(True)

        def _failed(self, message: str):
            self.status.setText(f"Could not compute: {message}")

        def _done(self, r, figs):
            text = r.report()
            self.report.setPlainText(text)
            self._fit_report_width(text.splitlines())
            previous = self.plot_tabs.currentIndex()
            while self.plot_tabs.count():
                w = self.plot_tabs.widget(0)
                self.plot_tabs.removeTab(0)
                w.deleteLater()
            for name, fig in figs.items():
                box = QtWidgets.QWidget()
                lay = QtWidgets.QVBoxLayout(box)
                lay.setContentsMargins(0, 0, 0, 0)
                canvas = FigureCanvasQTAgg(fig)
                lay.addWidget(NavigationToolbar2QT(canvas, box))
                lay.addWidget(canvas, 1)
                self.plot_tabs.addTab(box, PLOT_TITLES[name])
            if 0 <= previous < self.plot_tabs.count():
                self.plot_tabs.setCurrentIndex(previous)
            self.figures = figs
            self.save_button.setEnabled(True)
            self.status.setText(f"{r.individual}: risk {r.risk:.4g} · population "
                                f"{r.population_risk:.4g} · λ = {r.relative_risk:.2f}")

        def save(self):
            folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Plot folder")
            if not folder:
                return
            fmt, ok = QtWidgets.QInputDialog.getItem(
                self, "Format", "File format:", ["png", "svg", "pdf"], 0, False)
            if not ok:
                return
            try:
                for name, fig in self.figures.items():
                    fig.savefig(os.path.join(folder, f"{name}.{fmt}"), dpi=150,
                                bbox_inches="tight")
            except OSError as e:
                QtWidgets.QMessageBox.warning(self, "Could not save", str(e))
                return
            self.status.setText(f"{len(self.figures)} plots saved: {folder}")

    return MFPanel()
