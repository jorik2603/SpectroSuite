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

## Dependencies

- NumPy
- Matplotlib
- SciPy
- cmcrameri
