# Bioinformatics Project

This project is a small Python desktop application for exploring neural gene data.
It fetches selected gene sequences from NCBI, searches for user-provided DNA motifs,
and simulates an EEG-like signal based on motif frequency.

## Main Files

- `app.py`: Tkinter application for gene lookup, motif analysis, and EEG simulation
- `Bioinformatics_Report.pdf`: full project report with methodology, implementation details, and discussion
- `Bioinformatics_Research.pdf`: supporting research material

## Requirements

- Python 3
- `biopython`
- `numpy`
- `matplotlib`

## Run

```bash
pip install biopython numpy matplotlib
python app.py
```

An internet connection is required because the application fetches sequence data from NCBI.
For the complete explanation of the project, refer to `Bioinformatics_Report.pdf`.
