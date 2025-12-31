# Installation Instructions

## Method 1: Install from GitHub (Recommended)

```bash
# Install directly from GitHub
pip install git+https://github.com/jrneilson/molecularsymmetry.git

# Or clone and install in development mode
git clone https://github.com/jrneilson/molecularsymmetry.git
cd molecularsymmetry
pip install -e .
```

## Method 2: Local installation

```bash
# Download/clone the repository, then:
cd molecularsymmetry
pip install .

# For development with optional dependencies:
pip install -e ".[dev,jupyter]"
```

## Method 3: Copy module (legacy)

Simply copy `molecular_symmetry.py` to your project directory or Python path.

```bash
# Only external dependency is NumPy
pip install numpy
```

## Using in conda environments

The package can be installed in any conda environment:

```bash
# Activate your environment
conda activate myenv

# Install the package
pip install git+https://github.com/jrneilson/molecularsymmetry.git

# Or install with optional dependencies
pip install "git+https://github.com/jrneilson/molecularsymmetry.git[jupyter]"
```

## Usage after installation

```python
from molecular_symmetry import get_point_group

# Test the installation
oh = get_point_group('Oh')
print("Installation successful!")
oh.print_character_table()
```