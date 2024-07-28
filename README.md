[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/I2-kIT4b) \# Instructions

Basic structure to create a report using Quarto, based on Pandoc's markdown and LaTeX.

Check the doc at <https://quarto.org/> and simply install Quarto.

The `docs` folder can be directly used by Github Pages. Example:

<https://julien-vitay.net/quarto-report/>

## Instructions

-   Html and pdf documents are in the "docs" folder;

-   Since we used the Bestsubset package, we decided to save the workspace in order to make the code reproducible. The main code files are "sim.R" and "nonzero.R" ("workspace.RData"); "dfplot.R" ("workspace_df.RData"). Both the workspace and the code can be found in the folder "code";

-   In the folders "fig5" and "fig5_t" are stored the output of the simulations (.RDS files) necessary to reproduce some of the plots

-   The plots are in the folder "plots" in "dcs";

-   Check "renv" folder for the virtual environment

## Attributions:

This template is from <https://github.com/vitay/quarto-report>

## TODO before submitting

-   compile the report and check the [compiled pdf](docs/report.pdf)
-   Rewrite this README with instructions on how to re-run your code. Your work should be **reproducible**. Check the virtual environment page on the website for more information.
