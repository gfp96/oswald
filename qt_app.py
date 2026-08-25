"""PySide6/pyqtgraph frontend for Oswald.

This module is intentionally separate from ``main.py``.  The original
Tkinter interface remains available while this frontend provides a faster,
modern plotting path for evaluation and gradual migration.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from .Functions import Read_BE_file, Read_terratek_signal, Update_df
from .signal_interp import (
    Find_start,
    Find_start_burst,
    Get_Max_AIC_velocity,
    Get_STALTA_AIC_velocity,
    Get_filtered_signal,
)


@dataclass
class LoadedSignals:
    """Signals and their time vectors aligned with one dataframe selection."""

    times: list[np.ndarray]
    inputs: list[np.ndarray]
    outputs: list[np.ndarray]


class AnalysisWorker(QtCore.QObject):
    """Run disk I/O and numerical analysis away from the Qt event loop."""

    # Keep all aligned arrays in one payload.  Qt validates this declaration
    # at runtime, so every emit must provide exactly one argument.
    loaded = QtCore.Signal(object)
    analysed = QtCore.Signal(object, object, object, object)
    failed = QtCore.Signal(str)
    finished = QtCore.Signal()

    def __init__(self, folder, selection, data_format, operation, analysis_args=None, encap_files=None):
        super().__init__()
        self.folder = folder
        self.selection = selection
        self.data_format = data_format
        self.operation = operation
        self.analysis_args = analysis_args
        self.encap_files = encap_files or {}

    @QtCore.Slot()
    def run(self):
        """Dispatch work from the QThread's event loop, never from the GUI."""
        # QThread.started invokes this method in the worker's thread.  The
        # worker emits results; it never modifies Qt widgets directly.
        if self.operation == "load":
            self.load()
        else:
            self.analyse(*self.analysis_args)

    @QtCore.Slot()
    def load(self):
        """Read selected wave files and return one payload to the GUI thread.

        The selection dataframe is already owned by ``MainWindow``.  Only the
        newly loaded arrays cross the thread boundary, packaged in
        ``LoadedSignals`` so the signal declaration and emission stay aligned.
        """
        try:
            times, inputs, outputs = [], [], []
            if self.data_format == "Navier_BE":
                # Navier files contain both the excitation (input) and the
                # receiving (rough) waveform used by the interpretation.
                for filename in self.selection.filename:
                    time, input_signal, output_signal = Read_BE_file(
                        self.folder, filename
                    )
                    times.append(np.asarray(time, dtype=float))
                    inputs.append(np.asarray(input_signal, dtype=float))
                    outputs.append(np.asarray(output_signal, dtype=float))
            else:
                # Terratek files provide the receiving waveform.  Their
                # start index is derived later from the encap-to-encap time.
                for filename in self.selection.filename:
                    time, output_signal = Read_terratek_signal(
                        os.path.join(self.folder, filename)
                    )
                    times.append(np.asarray(time, dtype=float))
                    inputs.append(np.array([], dtype=float))
                    outputs.append(np.asarray(output_signal, dtype=float))
            # `loaded` is Signal(object), therefore this must be one argument.
            self.loaded.emit(LoadedSignals(times, inputs, outputs))
        except Exception as error:
            self.failed.emit(f"Could not load signals: {error}")
        finally:
            self.finished.emit()

    @QtCore.Slot(str, float, object, object, object, object)
    def analyse(self, method, length, signals, freq_range_p, freq_range_s, periods):
        """Run the selected interpretation method in the worker thread."""
        try:
            times = signals.times
            rough = signals.outputs
            starts = np.zeros(len(rough), dtype=int)
            for index, signal in enumerate(signals.inputs):
                # Estimate the beginning of the useful receiving signal before
                # filtering or AIC.  This is the same distinction as in main.py:
                # Terratek uses a reference travel time, bursts use their period
                # count, and ordinary BE signals use extrema in the input trace.
                if self.data_format == "Terratek":
                    dt = times[index][1] - times[index][0]
                    travel_time = 4.6e-6 if self.selection.isvp.iloc[index] else 7.1e-6
                    starts[index] = int(round(travel_time / dt))
                elif self.selection.Burst.iloc[index]:
                    starts[index] = Find_start_burst(signal, periods)
                else:
                    starts[index] = Find_start(signal)

            # Filter/Show CC only produce cleaned traces.  The AIC methods also
            # return arrival indices and velocities for the velocity panel.
            if method == "Filter":
                result = Get_filtered_signal(
                    times, starts, rough, self.selection.isvp,
                    self.selection.infreq, freq_range_p, freq_range_s,
                )
                self.analysed.emit(None, None, result, starts)
            elif method == "Max_AIC":
                result = Get_Max_AIC_velocity(
                    times, starts, rough, length, self.selection.infreq, freq_range_p
                )
                self.analysed.emit(result[0], result[1], result[2], starts)
            elif method == "STA/LTA_AIC":
                result = Get_STALTA_AIC_velocity(
                    times, starts, rough, length, self.selection.infreq, freq_range_p
                )
                self.analysed.emit(result[0], result[1], result[3], starts)
            else:
                self.analysed.emit(None, None, rough, starts)
        except Exception as error:
            self.failed.emit(f"Analysis failed: {error}")
        finally:
            self.finished.emit()


class MainWindow(QtWidgets.QMainWindow):
    """Interactive PySide6 replacement for the three Matplotlib canvases."""

    def __init__(self, folder, datafile):
        super().__init__()
        self.folder = Path(folder)
        self.datafile = Path(datafile)
        # `data` is the master table.  `selection` is only a filtered view;
        # edits are always written back to `data` before saving.
        self.signals = None
        self.selection = None
        self.filtered = None
        self.worker_thread = None
        self.worker = None
        self.setWindowTitle("Oswald | Wave propagation analysis")
        self.resize(1500, 950)
        self._load_database()
        self._build_ui()
        self._refresh_filter_values()
        self._refresh_selection()

    def _load_database(self):
        """Normalize the experimental table once, before any plotting."""
        # Spreadsheet input is intentionally normalized at the boundary so
        # the rest of the application can rely on stable dtypes and names.
        data = pd.read_excel(self.datafile)
        data = data.replace([np.inf, -np.inf], np.nan)
        data = data[(data.Sweep == False) & (data.freqlev >= 0)]
        self.data = data.drop_duplicates(
            subset={"isvp", "freqlev", "stageno"}, keep="first"
        ).copy()
        for column in ("stageno", "pelev", "pflev"):
            self.data[column] = self.data[column].astype(int)
        self.data["freqlev"] = self.data["freqlev"].astype(float)
        for column in ("isvp", "Sweep", "Burst", "valid", "valid_aicmax", "valid_SLA"):
            self.data[column] = self.data[column].astype(bool)
        # Keep all user and automatic results in the master table so filtering
        # and saving do not lose values from another stage.
        for column, default in (("Grade", -1), ("vel_manual", -1.0),
                                ("arrival_manual", np.nan),
                                ("vel_aicmax", np.nan), ("arrival_aicmax", np.nan),
                                ("vel_SLA", np.nan), ("arrival_SLA", np.nan)):
            Update_df(self.data, column, default)

    def _build_ui(self):
        """Create a compact control rail and three scalable plot panels."""
        # The grid gives the overview and velocity plot equal top-level
        # anchors while the detail view spans the full width below them.
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QGridLayout(central)
        controls = QtWidgets.QFormLayout()
        layout.addLayout(controls, 0, 0, 3, 1)

        self.pe = self._combo(["Any"])
        self.pf = self._combo(["Any"])
        self.stage = self._combo(["Any"])
        self.pe.currentTextChanged.connect(self._filters_changed)
        self.pf.currentTextChanged.connect(self._filters_changed)
        self.stage.currentTextChanged.connect(self._refresh_selection)
        controls.addRow("p' [kPa]", self.pe)
        controls.addRow("pf [kPa]", self.pf)
        controls.addRow("Stage", self.stage)

        self.format = QtWidgets.QComboBox()
        self.format.addItems(["Navier_BE", "Terratek"])
        controls.addRow("File format", self.format)
        self.method = QtWidgets.QComboBox()
        self.method.addItems(["None", "Filter", "Max_AIC", "STA/LTA_AIC", "Show CC"])
        controls.addRow("Method", self.method)
        self.method.currentTextChanged.connect(self._method_changed)
        self.frequency = self._combo([])
        self.frequency.setEnabled(False)
        self.frequency.currentIndexChanged.connect(self._frequency_changed)
        controls.addRow("Frequency [kHz]", self.frequency)
        self.length = QtWidgets.QDoubleSpinBox()
        self.length.setRange(1e-6, 10.0)
        self.length.setValue(0.193)
        self.length.setSuffix(" m")
        controls.addRow("BE length", self.length)
        self.periods = QtWidgets.QSpinBox()
        self.periods.setRange(1, 1000)
        self.periods.setValue(21)
        # This control only affects burst records.  It is updated whenever the
        # current stage selection changes, rather than forcing users to infer
        # whether the value is relevant.
        self.periods.setEnabled(False)
        controls.addRow("Burst periods", self.periods)

        self.load_button = QtWidgets.QPushButton("Load selected signals")
        self.load_button.clicked.connect(self._load_signals)
        controls.addRow(self.load_button)
        self.analyse_button = QtWidgets.QPushButton("Analyse")
        self.analyse_button.clicked.connect(self._analyse)
        self.analyse_button.setEnabled(False)
        controls.addRow(self.analyse_button)
        self.save_button = QtWidgets.QPushButton("Save results")
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self._save)
        controls.addRow(self.save_button)
        self.plot_button = QtWidgets.QPushButton("Plot velocity results")
        self.plot_button.setEnabled(False)
        self.plot_button.clicked.connect(self._draw_velocity)
        controls.addRow(self.plot_button)
        self.reset_button = QtWidgets.QPushButton("Reset pressure filters")
        self.reset_button.clicked.connect(self._reset_filters)
        controls.addRow(self.reset_button)
        self.status = QtWidgets.QLabel("Choose a stage and load its signals")
        self.status.setWordWrap(True)
        controls.addRow(self.status)

        pg.setConfigOptions(background="#101820", foreground="#e8edf2", antialias=True)
        self.main_plot = pg.PlotWidget(title="Stage waveform overview")
        self.velocity_plot = pg.PlotWidget(title="Velocity by frequency")
        self.detail = pg.GraphicsLayoutWidget()
        self.p_plot = self.detail.addPlot(row=0, col=0, title="Compression wave (left click)")
        self.s_plot = self.detail.addPlot(row=1, col=0, title="Shear wave (right click)")
        for plot in (self.main_plot, self.velocity_plot, self.p_plot, self.s_plot):
            plot.showGrid(x=True, y=True, alpha=0.2)
        self.main_plot.invertY(True)
        layout.addWidget(self.main_plot, 0, 1)
        layout.addWidget(self.velocity_plot, 0, 2)
        layout.addWidget(self.detail, 1, 1, 2, 2)
        self.main_plot.scene().sigMouseClicked.connect(self._overview_clicked)
        self.detail.scene().sigMouseClicked.connect(self._detail_clicked)

        # A normal Tab press is consumed by child widgets such as combo boxes.
        # This shortcut provides a dependable next-frequency command while
        # preserving Tab behavior when focus is in the main window itself.
        self.next_shortcut = QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Right"), self)
        self.next_shortcut.activated.connect(self._next_frequency)
        self.grade_shortcuts = []
        for grade in range(6):
            shortcut = QtGui.QShortcut(QtGui.QKeySequence(str(grade)), self)
            shortcut.setContext(QtCore.Qt.ShortcutContext.WindowShortcut)
            shortcut.activated.connect(lambda grade=grade: self._grade(grade))
            self.grade_shortcuts.append(shortcut)

    @staticmethod
    def _combo(values):
        """Create a combo box from display values."""
        combo = QtWidgets.QComboBox()
        combo.addItems(values)
        return combo

    def _refresh_selection(self):
        """Apply both pressure filters and the selected stage."""
        # Pressure choices constrain the stage choice; the stage then applies
        # one additional filter to the already reduced metadata table.
        result = self.data
        if self.pe.currentText() != "Any":
            result = result[result.pelev == int(self.pe.currentText())]
        if self.pf.currentText() != "Any":
            result = result[result.pflev == int(self.pf.currentText())]
        stage = self.stage.currentText()
        if stage != "Any":
            result = result[result.stageno == int(stage)]
        self.selection = result.copy()
        self.periods.setEnabled(bool(self.selection.get("Burst", pd.Series(dtype=bool)).any()))
        self.status.setText(f"{len(self.selection)} recordings selected")

    def _refresh_filter_values(self):
        """Populate pressure and stage selectors without triggering recursion."""
        self.pe.addItems(str(value) for value in sorted(self.data.pelev.unique()))
        self.pf.addItems(str(value) for value in sorted(self.data.pflev.unique()))
        self._update_stage_values()

    def _update_stage_values(self):
        """Limit stages to those compatible with the pressure selections."""
        # Block signals while rebuilding the combo to avoid recursive updates
        # when a previously selected stage is no longer valid.
        result = self.data
        if self.pe.currentText() != "Any":
            result = result[result.pelev == int(self.pe.currentText())]
        if self.pf.currentText() != "Any":
            result = result[result.pflev == int(self.pf.currentText())]
        current = self.stage.currentText()
        values = ["Any"] + [str(value) for value in sorted(result.stageno.unique())]
        self.stage.blockSignals(True)
        self.stage.clear()
        self.stage.addItems(values)
        self.stage.setCurrentText(current if current in values else "Any")
        self.stage.blockSignals(False)

    def _filters_changed(self):
        """Refresh stages and selected rows after a pressure change."""
        self._update_stage_values()
        self._refresh_selection()

    def _reset_filters(self):
        """Restore the complete dataset and all available stages."""
        self.pe.setCurrentText("Any")
        self.pf.setCurrentText("Any")
        self._update_stage_values()
        self.stage.setCurrentText("Any")
        self._refresh_selection()

    def _method_changed(self):
        """Re-run the selected display/interpretation method when loaded."""
        if self.signals is not None:
            self._analyse()

    def _load_signals(self):
        """Start asynchronous loading so the window remains responsive."""
        self._start_worker("load")

    def _analyse(self):
        """Start asynchronous interpretation for the loaded selection."""
        if self.signals is None:
            return
        self._start_worker("analyse")

    def _start_worker(self, operation):
        self.load_button.setEnabled(False)
        self.analyse_button.setEnabled(False)
        self.worker_thread = QtCore.QThread(self)
        analysis_args = None
        if operation == "analyse":
            analysis_args = (
                self.method.currentText(), self.length.value(), self.signals,
                np.array([100.0, 200_000.0]), np.array([100.0, 100_000.0]),
                self.periods.value(),
            )
        self.worker = AnalysisWorker(
            self.folder, self.selection, self.format.currentText(), operation, analysis_args
        )
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.loaded.connect(self._signals_loaded)
        self.worker.analysed.connect(self._analysis_complete)
        self.worker.failed.connect(self._show_error)
        self.worker.finished.connect(self._worker_finished)
        self.worker_thread.finished.connect(self.worker.deleteLater)
        self.worker_thread.start()

    @QtCore.Slot(object)
    def _signals_loaded(self, signals):
        self.signals = signals
        self.filtered = None
        self.analysis = None
        self.analysis_curves = None
        self.starts = None
        self.current_index = 0
        self.frequency.blockSignals(True)
        self.frequency.clear()
        self.frequency.addItems(str(value) for value in sorted(self.selection.freqlev.unique()))
        self.frequency.setEnabled(True)
        self.frequency.blockSignals(False)
        self._draw_overview()
        self._draw_detail(0)
        self._draw_velocity()
        self.plot_button.setEnabled(True)
        self.analyse_button.setEnabled(True)
        self.status.setText(f"Loaded {len(signals.outputs)} recordings")

    @QtCore.Slot(object, object, object, object)
    def _analysis_complete(self, result, curves, filtered, starts):
        """Store results and refresh every visualization from the dataframe."""
        self.analysis = result
        self.analysis_curves = curves
        self.filtered = filtered
        self.starts = starts
        self._store_analysis(result)
        self._draw_overview()
        self._draw_velocity()
        self._draw_detail(self.current_index)
        self.save_button.setEnabled(True)
        self.plot_button.setEnabled(True)
        self.analyse_button.setEnabled(True)
        self.status.setText("Analysis complete")

    def _store_analysis(self, result):
        """Copy automatic arrivals and velocities into the master dataframe."""
        if result is None:
            return
        method = self.method.currentText()
        for row_number, (_, row) in enumerate(self.selection.iterrows()):
            filename = row.filename
            if method == "Max_AIC":
                values = (result.iloc[row_number].vel, result.iloc[row_number].a_time)
                self.data.loc[self.data.filename == filename, ["vel_aicmax", "arrival_aicmax"]] = values
            elif method == "STA/LTA_AIC":
                velocity = result.iloc[row_number].vel1 if row.isvp else result.iloc[row_number].vel2
                arrival = result.iloc[row_number].a1_time if row.isvp else result.iloc[row_number].a2_time
                self.data.loc[self.data.filename == filename, ["vel_SLA", "arrival_SLA"]] = velocity, arrival
        self.selection = self.data.loc[self.selection.index].copy()

    def _draw_overview(self):
        """Spread channels along X and retain high-resolution arrays for zooming."""
        # A frequency value is not a useful X coordinate for a waveform stack:
        # P and S records may share it.  Each row therefore receives a fixed
        # channel offset, with amplitude normalized around that offset.
        self.main_plot.clear()
        if self.signals is None or self.selection.empty:
            return
        signals = self.filtered if self.filtered is not None else self.signals.outputs
        traces = [signal for signal in signals if signal.size]
        if not traces:
            return
        amplitude = max(
            np.nanmax(np.abs(signal)) for signal in traces
        ) or 1.0
        self.overview_centres = []
        self.overview_labels = []
        labelled_frequencies = set()
        for index, (time, signal, isvp) in enumerate(zip(self.signals.times, signals, self.selection.isvp)):
            if not signal.size:
                continue
            centre = index * 2.4
            self.overview_centres.append(centre)
            # TextItem uses plot coordinates, so the label follows the trace
            # when the user zooms or pans instead of being painted on the UI.
            frequency = float(self.selection.freqlev.iloc[index])
            if frequency not in labelled_frequencies:
                label = pg.TextItem(str(frequency), color="#e8edf2", anchor=(0.5, 1.0))
                label.setPos(centre, float(np.nanmax(time) * 1000))
                self.main_plot.addItem(label)
                self.overview_labels.append(label)
                labelled_frequencies.add(frequency)
            pen = pg.mkPen("#55c2ff" if isvp else "#f2b134", width=1.2)
            # A dashed excitation trace plus a solid receiving trace matches
            # the two traces shown in the original top-left figure.
            input_signal = self.signals.inputs[index]
            if input_signal.size:
                self.main_plot.plot(input_signal / amplitude + centre, time[:input_signal.size] * 1000,
                                    pen=pg.mkPen("#8b949e", style=QtCore.Qt.PenStyle.DashLine))
            curve = self.main_plot.plot(signal / amplitude + centre, time * 1000, pen=pen)
            # Keep the original arrays for detail plots, but let pyqtgraph
            # draw only visible extrema when many samples occupy one pixel.
            curve.setClipToView(True)
            curve.setDownsampling(auto=True, method="peak")
        self.main_plot.setLabel("left", "Time", units="ms")
        self.main_plot.setLabel("bottom", "Recording channel; P blue, S amber")

    def _draw_velocity(self):
        """Show manual and automatic P/S velocities against frequency."""
        self.velocity_plot.clear()
        if self.selection.empty:
            return
        styles = {"vel_manual": ("#ec6a5e", "Manual", "o"),
                  "vel_aicmax": ("#55c2ff", "Max AIC", "t"),
                  "vel_SLA": ("#f2b134", "STA/LTA AIC", "x"),
                  "vel_CC": ("#a78bfa", "Cross-correlation", "s")}
        for column, (color, label, symbol) in styles.items():
            if column not in self.selection:
                continue
            for isvp, suffix in ((True, " P"), (False, " S")):
                rows = self.selection[(self.selection.isvp == isvp) & (self.selection[column] > 0)]
                if not rows.empty:
                    self.velocity_plot.plot(rows.freqlev.to_numpy(), rows[column].to_numpy(),
                                            pen=None, symbol=symbol, symbolBrush=color,
                                            symbolPen=color, symbolSize=9, name=label + suffix)
        self.velocity_plot.setLabel("left", "Velocity", units="m/s")
        self.velocity_plot.setLabel("bottom", "Frequency", units="kHz")
        self.velocity_plot.setYRange(0, 2000, padding=0)

    def _overview_clicked(self, event):
        """Select the nearest offset channel in the overview."""
        if self.signals is None or event.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        position = self.main_plot.plotItem.vb.mapSceneToView(event.scenePos())
        centres = np.asarray(getattr(self, "overview_centres", []))
        if centres.size:
            row_index = int(np.argmin(np.abs(centres - position.x())))
            frequency = str(self.selection.freqlev.iloc[row_index])
            self.frequency.setCurrentText(frequency)

    def _frequency_changed(self, index):
        """Render the frequency selected in the combo box."""
        if self.signals is not None and 0 <= index < self.frequency.count():
            self._draw_detail(index)

    def _next_frequency(self):
        """Advance to the next unique frequency without moving widget focus."""
        if self.signals is not None and self.frequency.count():
            self.frequency.setCurrentIndex((self.current_index + 1) % self.frequency.count())

    def _draw_detail(self, index):
        """Draw separate P and S traces with start and saved-arrival markers."""
        if self.signals is None or self.selection.empty:
            return
        index = max(0, min(index, self.frequency.count() - 1))
        frequency = float(self.frequency.itemText(index))
        self.current_index = index
        self.p_plot.clear()
        self.s_plot.clear()
        self.detail_row_indices = {}
        self.detail_views = getattr(self, "detail_views", {})
        for plot, isvp in ((self.p_plot, True), (self.s_plot, False)):
            # One frequency can have two rows, one P and one S.  Locate each
            # independently so both subplots compare the same frequency.
            matches = self.selection[(self.selection.freqlev == frequency) & (self.selection.isvp == isvp)]
            if matches.empty:
                plot.setTitle(f"No {'compression' if isvp else 'shear'} wave at this frequency")
                continue
            row_index = self.selection.index.get_loc(matches.index[0])
            self.detail_row_indices[isvp] = row_index
            row = matches.iloc[0]
            time = self.signals.times[row_index] * 1000
            signal = self.signals.outputs[row_index]
            if self.filtered is not None:
                signal = self.filtered[row_index]
            input_signal = self.signals.inputs[row_index]
            if input_signal.size:
                # The input and received signals can differ by orders of
                # magnitude.  Overlay the input in a linked ViewBox so its
                # right-hand scale is independent of the output scale.
                input_view = self.detail_views.get(plot)
                if input_view is None:
                    input_view = pg.ViewBox()
                    self.detail_views[plot] = input_view
                    plot.scene().addItem(input_view)
                    plot.showAxis("right")
                    plot.getAxis("right").linkToView(input_view)
                    input_view.setXLink(plot)
                    plot.vb.sigResized.connect(
                        lambda view=input_view, source=plot: view.setGeometry(source.vb.sceneBoundingRect())
                    )
                input_view.clear()
                input_view.addItem(pg.PlotDataItem(
                    time[:input_signal.size], input_signal,
                    pen=pg.mkPen("#8b949e", style=QtCore.Qt.PenStyle.DashLine),
                ))
                input_view.setYRange(float(np.nanmin(input_signal)), float(np.nanmax(input_signal)), padding=0.1)
                plot.getAxis("right").setLabel("Input amplitude", color="#8b949e")
            plot.plot(time, signal, pen=pg.mkPen("#55c2ff" if isvp else "#f2b134", width=1.4))
            if self.starts is not None and self.starts[row_index] < len(time):
                plot.addLine(x=time[self.starts[row_index]],
                             pen=pg.mkPen("#e8edf2", style=QtCore.Qt.PenStyle.DashLine))
            for column, color in (("arrival_manual", "#ec6a5e"),
                                  ("arrival_aicmax", "#55c2ff"),
                                  ("arrival_SLA", "#f2b134")):
                if pd.notna(row[column]):
                    plot.addLine(x=float(row[column]) * 1000, pen=pg.mkPen(color, width=2))
            plot.setLabel("left", "Amplitude")
            plot.setLabel("bottom", "Time", units="ms")
            plot.setTitle(f"{'P' if isvp else 'S'} wave | {row.freqlev:g} kHz | click to pick")

    def _detail_clicked(self, event):
        """Store a left-click P or right-click S arrival and its velocity."""
        if self.signals is None or self.selection.empty:
            return
        scene_pos = event.scenePos()
        plot = self.p_plot if self.p_plot.sceneBoundingRect().contains(scene_pos) else self.s_plot
        isvp = plot is self.p_plot
        expected_button = QtCore.Qt.MouseButton.LeftButton if isvp else QtCore.Qt.MouseButton.RightButton
        if event.button() != expected_button:
            return
        row_index = self.detail_row_indices.get(isvp)
        if row_index is None:
            return
        row = self.selection.iloc[row_index]
        if bool(row.isvp) != isvp:
            return
        # pyqtgraph reports the click in milliseconds because the detail plots
        # use milliseconds on X; the stored arrival remains in seconds to
        # match the signal-processing functions.
        arrival = plot.vb.mapSceneToView(scene_pos).x() / 1000
        start = self.signals.times[row_index][self.starts[row_index]]
        if arrival <= start:
            self.status.setText("Arrival must be later than the detected start")
            return
        # Start-to-start velocity is the known BE tip-to-tip distance divided
        # by the manually picked propagation time.
        velocity = self.length.value() / (arrival - start)
        filename = row.filename
        self.data.loc[self.data.filename == filename, ["vel_manual", "arrival_manual"]] = velocity, arrival
        self.selection = self.data.loc[self.selection.index].copy()
        self._draw_velocity()
        self._draw_detail(self.current_index)
        self.save_button.setEnabled(True)
        self.status.setText(f"Manual {'P' if isvp else 'S'} pick: {velocity:.1f} m/s")

    def keyPressEvent(self, event):
        """Advance frequency with Tab and grade the current signal with 0-5."""
        if event.key() == QtCore.Qt.Key.Key_Tab and self.signals is not None:
            self.frequency.setCurrentIndex((self.current_index + 1) % self.frequency.count())
            event.accept()
            return
        if QtCore.Qt.Key.Key_0 <= event.key() <= QtCore.Qt.Key.Key_5:
            grade = event.key() - QtCore.Qt.Key.Key_0
            filename = self.selection.iloc[self.current_index].filename
            self.data.loc[self.data.filename == filename, "Grade"] = grade
            if grade == 0:
                self.data.loc[self.data.filename == filename, "valid"] = False
            self.selection = self.data.loc[self.selection.index].copy()
            self.status.setText(f"Grade {grade} stored")
            event.accept()
            return
        super().keyPressEvent(event)

    def _save(self):
        """Write all current manual and automatic results back to Excel."""
        try:
            self.data.to_excel(self.datafile, index=False)
            self.status.setText(f"Saved results to {self.datafile.name}")
        except Exception as error:
            self._show_error(f"Could not save results: {error}")

    def _worker_finished(self):
        self.load_button.setEnabled(True)
        if self.signals is not None:
            self.analyse_button.setEnabled(True)
        if self.worker_thread is not None:
            self.worker_thread.quit()
            self.worker_thread.wait()

    def _show_error(self, message):
        self.status.setText(message)
        QtWidgets.QMessageBox.critical(self, "Oswald", message)


def main():
    """Launch the Qt frontend without creating a window during import."""
    app = QtWidgets.QApplication(sys.argv)
    folder = QtWidgets.QFileDialog.getExistingDirectory(None, "Select signal folder")
    if not folder:
        return 0
    datafile, _ = QtWidgets.QFileDialog.getOpenFileName(
        None, "Select experimental database", folder, "Excel files (*.xlsx);;All files (*)"
    )
    if not datafile:
        return 0
    window = MainWindow(folder, datafile)
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())