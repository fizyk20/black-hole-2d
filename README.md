# Black Hole 2D

## What is it?

It's a simple 2-dimensional simulator of flight near a black hole.

The ship moves according to laws of General Relativity.

## Controls

You can pan the view by dragging it with the mouse. Right-clicking and dragging zooms the view.

Keyboard controls:

- W/S/A/D: burn engines forward/backwards/left/right
- U/O: rotate the ship left/right
- R/T: decrease/increase time warp
- Shift: holding it makes the thrust accumulate, that is - increase while the thrust controls are held

Buttons "Start", "Stop" and "Reset" control the simulation.

The parameters of the black hole, as well as the position of the ship, can be adjusted using the text
boxes on the right.

By default, the simulation is run according to time at infinity - in such a setting it is impossible
to fall into the black hole, as it would take infinitely long to do so.
"Use on-board time for simulation" makes the simulation run in on-board time, which can differ a lot
from time at infinity when moving at high speeds or being close to/in the black hole.

## Compilation

Requires Qt 5. To compile:

```
$ qmake
$ make
```
