# pyfd

![pyfd logo](/img/logo.png)

A simple demonstration of finite differences for plotting potentials arising from 2D charge distributions, written for my EM3 class.

### Features:
not very many! Uses fixed potential boundaries.
Written in python.

Indexing is facaded behind compass directions for ease of writing and writing (avoids code with +/-1 all through it)

### Usage:

Download this repository, then run from the command line (or from your IDE)

```python fd.py```

Assumes you have python3 installed, along with numpy and matplotlib (which you will if you use a complete distribution like Anaconda, and if not, can be easily rectified e.g. using ```pip install numpy``` at the command line.

### Going further

#### Geometries

Implement your own geometry construction functions either by extending the existing Grid class, or by adding a new class that calls Grid.fixV().

Suggestions: circles, squares, rectangles, triangles, letters of the alphabet, cartoons

#### Accuracy

You might like to look more closely at convergece behaviour, and compare results to those obtained by hand, to understand just how much error was introduced by our assumption of an ideal capacitor. 

#### Validation
I've not done this yet, because things "look about right" - at least, so far! I've left it as an exercise for you. That might seem a bit cheeky, but actually in industry and research, any time you get a new simulation tool, you undertake a validation exercise to compare it to measurements (directly or indirectly).


