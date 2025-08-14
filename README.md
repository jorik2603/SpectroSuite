# SpectroSuite

SpectroSuite is an interactive Python GUI for visualizing and comparing 2D spectroscopy datasets.

## Features

- **Multiple Analysis Modes**: Choose between viewing a single dataset, comparing two, or comparing up to four.
- **Interactive Slicing**: Use sliders or type exact axis values to view 1D slices.
- **Data Normalization**: Normalize data to its maximum value or scale to a `[0, 1]` or `[-1, 1]` range.
- **Advanced Plotting**: Includes data windowing, Savitzky-Golay smoothing, and scientific colormaps from `cmcrameri`.

## Installation

You can install SpectroSuite directly from the root directory using pip:

```bash
pip install .
```

This will install the package and all required dependencies.

## Usage

Once installed, you can launch the application from your terminal by running:

```bash
spectrosuite
```

This will open a dialog asking you to select an analysis mode.

## Supported Data Formats 📁

The application supports several data formats. For the program to correctly interpret the axes and data, your files should be structured as follows:

### Specific `.dat` Format (Recommended)

This is the primary supported format, as it includes the x and y axis values directly in the file.

- The **first line** must be a header containing the **x-axis** values (e.g., wavelength). It should start with a `0` or placeholder, followed by comma-separated numerical values.
- **Each subsequent line** should start with a single **y-axis** value (e.g., time), followed by the comma-separated data values for that row.

**Example:**
```
0,800,805,810,815     # Line 1: Placeholder, followed by x-axis values (wavelength)
0.1,10,12,15,14       # Line 2: y-axis value (time), followed by data
0.2,25,28,32,30       # Line 3: y-axis value (time), followed by data
0.3,41,45,50,48       # etc...
```

### General Text Files (`.txt`, `.csv`)

These files should contain **only the numerical data matrix**. The values should be separated by commas (or other whitespace recognized by `numpy.loadtxt`). The program will treat the row and column numbers as the y and x axes, respectively.

**Example:**
```
10,12,15,14
25,28,32,30
41,45,50,48
```

### NumPy Binary Files (`.npy`)

You can load data saved as a standard NumPy binary array file. The file should contain a single 2D NumPy array. As with `.txt` files, the axes will be inferred from the array indices.


## Dependencies

- NumPy
- Matplotlib
- SciPy
- cmcrameri
