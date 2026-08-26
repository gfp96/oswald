# oswald: OpenSource Wave Analysis for Laboratory soil Dynamics

## Overview
Oswald is a graphical interface designed for the interpretation large numbers of wave propagation measurements in porous media.

## Installation
Oswald is currently only available through cloning the github repository.
```ShellSession
cd C:/path/to/installation/folder
git clone https://github.com/gfp96/oswald
```

Then the GUI is launched by running ***main.py***

### PySide6/pyqtgraph preview

The newer visual layer can be launched from the workspace root after installing
the Qt dependencies:

```ShellSession
cd 03_lab_management
python -m pip install -r oswald/requirements-qt.txt
python -m oswald.qt_app
```

This frontend keeps the existing `main.py` Tkinter application available. It
loads wave files and runs interpretation in a worker thread so the interface
remains responsive, and uses pyqtgraph for interactive waveform rendering.

## How does it work
Oswald enables the user to manually pick arrival times using the Start-to-Start method. It also allows the use of two automated wave arrival picking methods: the MAIC and SLA methods described in [1]

Oswald requires two separate inputs:
- database of recorded signals in xlsx format
- 1 file per recording containing the rough data

### Input database
The input database must have a specific column names to ensure the filtering tools function correctly. Oswald is currently unable to manage duplicate entries in this database. Only the first entry will be kept.
Necessary data columns:
- 

### Keyboard bindings

The following keyboard bindings are available in the graphical interfaces:

| Binding | Action |
| --- | --- |
| `0` to `5` | Assigns the corresponding quality grade to the currently selected signal. In the Qt interface, the grade is kept locally until **Save results** is pressed. |
| `Tab` | Moves to the next frequency in the signal list in the Tkinter interface and, when focus permits, in the Qt interface. |
| `Ctrl` + `Right Arrow` | Moves to the next frequency in the Qt interface, including when focus is in another control. |
| Left mouse button on the P-wave detail plot | Selects the compression-wave arrival, calculates its Start-to-Start velocity, and stores the manual pick. |
| Left mouse button on the S-wave detail plot | Selects the shear-wave arrival, calculates its Start-to-Start velocity, and stores the manual pick. |

Manual arrivals and grades can be reviewed in the detail and velocity plots. Use
**Save results** to write the current results back to the experimental database.

### Code structure

The project is currently organized as a small source-tree application:

- `main.py` contains the original Tkinter GUI, including controls, Matplotlib
	figures, manual picking, grading, filtering, and saving.
- `qt_app.py` contains the PySide6/pyqtgraph interface. It provides the newer
	visual layer, background workers for file loading and analysis, interactive
	waveform plots, velocity plots, and manual picking.
- `signal_interp.py` contains the signal-processing algorithms: start
	detection, filtering, Max-AIC interpretation, STA/LTA-AIC interpretation,
	and numerical helper functions.
- `Functions.py` contains shared file readers, dataframe helpers, plotting
	defaults, and compatibility utilities used by the GUI layers.
- `tests/` contains automated tests for the signal-processing functions and
	edge cases.
- `requirements-qt.txt` lists the PySide6 and pyqtgraph dependencies required
	by the newer interface.
- `__init__.py` makes the source directory importable as the `oswald` package.

The intended long-term architecture is to keep `signal_interp.py` and the I/O
helpers independent of the GUI, allowing the Tkinter frontend to remain
available while the Qt frontend and package API mature.



## References
[1] Flood-Page, Guillaume, Luc Boutonnier, and Jean-Michel Pereira. 2024. ‘Application of the Akaike Information Criterion to the Interpretation of Bender Element Tests’. Soil Dynamics and Earthquake Engineering 177 (February): 108373. https://doi.org/10.1016/j.soildyn.2023.108373.

