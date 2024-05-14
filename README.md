# Guide for MERFISH analysis (from MERSCOPE data)
This is a guide to utilizing custom analysis scripts for data generated from MERSCOPE instruments. A workflow will be updated below to generate segmented data as well as integration of useful tools/scripts for data analysis.

## General workflow:

1. Generate MERSCOPE data
2. Use vizgen-postprocessing tool to resegment data (Cellpose)
3. Perform clustering analysis of identified cells
4. Reload segments in .vzg file and into MERSCOPE visualizer

## Useful links:

- https://github.com/Vizgen/vizgen-postprocessing
- https://vizgen.github.io/vizgen-postprocessing/index.html

