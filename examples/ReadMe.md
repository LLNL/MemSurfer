#### MemSurfer Examples

Two simple examples are included to demonstrate the functionality of MemSurfer.

* **ex_simple.py:** This script reads a point set that represents a 2D sinusoidal surface with added noise in z-dimension.
* **ex_bilayer.py:** This script demonstrates a more realistic use case of a bilayer membrane consisting of three types of lipids. To read and parse the input data, the following are used.
  - `MDAnalysis`. Please install as `pip install mdAnalysis`.
  - `MDReader`. We include an open-source utility for parsing the data. Credit is due to Manuel N Melo for creating MDReader. Please find the latest version [here](https://github.com/mnmelo/mdreader).
  - `lipiddefs`. We include a custom utility to parse the example data.

Both the example scripts create several output files.
* 

### License

MemSurfer is released under GPU-3.0 license. See the `LICENSE` file for details.

*`LLNL-CODE-763493`*
