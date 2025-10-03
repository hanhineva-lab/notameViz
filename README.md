# notameViz - Non-targeted metabolomics visualization

The notameViz package is the notame package that comes in handy when 
visualizing non-targeted metabolomics data analysis. Includes quality 
control visualizations, feature-wise visualizations and results visualizations.
Check out the [website](https://hanhineva-lab.github.io/notame/) for more on 
notame.

## Installation and getting started

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("notameViz")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("notameViz")
```

If you find any bugs or other things to fix, please submit an issue on GitHub! 
All contributions to notame are always welcome!