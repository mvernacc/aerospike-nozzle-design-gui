Aerospike Nozzle Design GUI
===========================

Graphical software for designing aerospike rocket nozzles.


Usage
-----
Run the software with the following terminal command:

```
python nozzle_gui.py
```

The GUI provides fields for entering the perameters of your propulsion system (e.g. chamber temperature, cahmber pressure, expansion ratio, etc.). After new numbers are entered in any field, the solver will be run. After the solver runs, the nozzle design results will be displayed in a text field on the righthand side of the GUI. A scond window will also pop up, which shows plots of spike radius, pressure, Isp, Mach number, and Temperature versus spike length.

You may save and load the inputs parameters with the "Save Design" and "Load Design" buttons.

You may export the nozzle geometry to a list of points in a text file with the "Export Spike+Shroud Shape to CAD" button.  The points text file can then be imported into SolidWorks CAD software. In SolidWorks, use the "Insert > Curve > Curve Through XYZ Points ..." command. Then revolve the curve around the y axis.


Installation
------------
You need a python environment with Tkinter, numpy, and matplotlib. On Linux, run in terminal:

```
sudo apt-get install python python-numpy python-scipy python-tk python-matplotlib
```

On Windows 7 I have tested [WinPython](https://winpython.github.io/), although other python distributions might work.


Solver
------
The nozzle design solver is based on an algorithm developed by C.C. Lee for NASA’s Marshal Space Flight Center. A description and derivation of the algorithm may be found in [Spike Contour Algorithm.docx](Spike Contour Algorithm.docx).
