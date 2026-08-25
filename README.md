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



## References
[1] Flood-Page, Guillaume, Luc Boutonnier, and Jean-Michel Pereira. 2024. ‘Application of the Akaike Information Criterion to the Interpretation of Bender Element Tests’. Soil Dynamics and Earthquake Engineering 177 (February): 108373. https://doi.org/10.1016/j.soildyn.2023.108373.

