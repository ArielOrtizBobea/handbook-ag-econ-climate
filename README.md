# Replication code for Ortiz-Bobea (2021), Handbook of Agricultural Economics

R code that reproduces the figures in

> Ortiz-Bobea, A. (2021). The empirical analysis of climate change impacts and adaptation in agriculture. In C. B. Barrett and D. R. Just (Eds.), *Handbook of Agricultural Economics*, Vol. 5, Chapter 76, pp. 3981–4073. Elsevier. https://doi.org/10.1016/bs.hesagr.2021.10.002

This is an updated version of the code archived with the chapter at the Cornell Institute for Social and Economic Research (CISER), https://doi.org/10.6077/fb1a-c376. The 2021 code relied on R packages that have since been removed from CRAN (rgdal, rgeos, maptools) or superseded. This version uses current packages and keeps the structure of the 2021 scripts. It downloads its data, so the repository holds only code.

## What changed

- sf and terra replace sp, rgdal, rgeos, maptools and raster for spatial data.
- fixest replaces lfe for fixed-effects regressions. Its small-sample corrections are set to match lfe, so the standard errors are the same as in the chapter.
- The world map in Figure 1 comes from Natural Earth (rnaturalearth), where the 2021 code used rworldmap.
- `0_download_data.R` downloads the data. A few files that providers no longer distribute, or that took a day of computing to produce, are attached to a release of this repository.
- Each figure sits in a marked block that can be re-run on its own (`make block name=fig8`).

The updated code reproduces the published figures, and its estimates and standard errors agree with those of the 2021 code to at least seven significant digits (see [How the update was checked](#how-the-update-was-checked)).

## Quick start

Requirements: R (the figures were produced with R 4.6.1), and `make` if you want to use the commands below. On Linux, sf and terra also need the GDAL, GEOS and PROJ system libraries (see https://r-spatial.github.io/sf/#installing).

```bash
git clone https://github.com/ArielOrtizBobea/handbook-ag-econ-climate.git
cd handbook-ag-econ-climate
make all
```

`make all` installs the R packages (see below), runs the six scripts in order and writes the figures to `figures/`. Alternatively, open `handbook_R/handbook_R.Rproj` in RStudio, run `renv::restore()` once, and run the scripts in the order they are numbered. The scripts expect `handbook_R` to be the working directory, which the RStudio project sets.

The default run downloads about 650 MB and takes about 10 minutes on a recent laptop, most of it in `1_weather_data.R`. It uses about 1.5 GB of disk space and up to 6 GB of memory.

### R packages

The version of every R package used to produce the figures is recorded in `handbook_R/renv.lock` with [renv](https://rstudio.github.io/renv/). `renv::restore()` (run by `make packages` and `make all`) installs these versions into a library inside the project, so your other R libraries are not affected. The packages come from the snapshot of CRAN taken on 1 October 2026 by Posit Package Manager, which provides pre-built packages for Windows, macOS and Linux, so the restore takes a few minutes. The lockfile was created with R 4.6.1. With another version of R, renv shows a warning and may have to build some packages from source, which requires compilers.

To use your own R library instead, delete `handbook_R/.Rprofile`. Each script then installs any package it is missing from CRAN. The main versions used are sf 1.1-3, terra 1.9-50, fixest 0.14.2, prism 0.3.0, spdep 1.4-2, splm 1.6-5 and data.table 1.18.6.1. Two version requirements matter: prism 0.3.0 or later, because PRISM changed its download service in 2025, and fixest 0.14 or later, for the names of the arguments of `ssc()`.

## Scripts

All scripts are in `handbook_R/`.

| Script | Content | Figures |
|---|---|---|
| `0_download_data.R` | Downloads the data | |
| `1_weather_data.R` | Land cover over the PRISM grid, the matrix that aggregates gridded data to counties, interpolation of weather station data, temperature exposure bins, aggregation to counties | 1–7 |
| `2_nonlinear_effects.R` | Effects of temperature on corn yields with step functions, natural cubic splines and Chebyshev polynomials | 8–10 |
| `3_time-varying_effects.R` | Effects that vary over the growing season, with a tensor-product spline | 11–12 |
| `4_spatial_dependence.R` | Standard errors under spatial dependence and a spatial error model | 13 |
| `5_robustness_checks.R` | 72 alternative specifications shown in a specification chart | 14 |
| `spec_chart_reproducible_example.R` | Examples of the specification chart function on the Harrison and Rubinfeld (1978) housing data | |
| `functions.R` | Functions used by the scripts | |

`_run_block.R` runs a single block (see below) and `_make_release_assets.R` documents how the release files were built from the CISER archive. Neither is needed to reproduce the figures.

## Data

| Data | Provider | Downloaded from |
|---|---|---|
| County boundaries, 2010, 1:20 million | US Census Bureau | census.gov |
| Country boundaries, 1:10 million | Natural Earth | naturalearthdata.com, with rnaturalearth |
| Daily and monthly gridded weather, August 2020 | PRISM Climate Group, Oregon State University | PRISM web service, with the prism package |
| Weather station list and records for August 2020 | NOAA, GHCN-Daily | release (files of May 2021) |
| County corn yields, 1981–2020 | USDA NASS Quick Stats | release (export of May 2021) |
| Share of each land cover class in each PRISM grid cell | computed from the NLCD 2016 (MRLC) | release |
| County-level PRISM weather, 1981–2020 | computed from PRISM | release |

The release files are attached to the release [`data-v1`](https://github.com/ArielOrtizBobea/handbook-ag-econ-climate/releases/tag/data-v1) of this repository. They are the files used for the chapter in 2021, copied from the CISER archive. The reasons for archiving each one:

- **County-level PRISM weather.** Producing these files with `1_weather_data.R` requires about 30,000 daily PRISM grids (about 120 GB of disk) and a day of computing. They can be rebuilt (see below).
- **Land cover shares.** They were computed from the NLCD 2016 land cover release of April 2019, which MRLC no longer distributes. Later NLCD releases revised the 2016 map.
- **Corn yields.** The Quick Stats API requires a key. The export can be repeated at https://quickstats.nass.usda.gov with Program SURVEY, Sector CROPS, Group FIELD CROPS, Commodity CORN, Category YIELD, Data Item "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE", Geographic Level COUNTY and Years 1981–2020, split into two periods because of the export size limit.
- **GHCN-Daily.** NOAA revises the archive continuously. In October 2026, the 2020 file had about 375 more stations reporting maximum temperature on August 16, 2020 than in May 2021, and revised values for 665 of the stations present in both. Setting `ghcn_current <- TRUE` in `0_download_data.R` (after deleting `data/GHCN`, if present) uses today's files, which changes Figures 1, 3 and 4 slightly.

Please cite the data providers when using these data. PRISM data: PRISM Climate Group, Oregon State University, https://prism.oregonstate.edu. GHCN-Daily: Menne, M. J., et al. (2012), https://doi.org/10.7289/V5D21VHZ.

## Rebuilding the county-level weather data

Set `full <- TRUE` at the top of `0_download_data.R` and of `1_weather_data.R`, then run both scripts. The first downloads daily minimum and maximum temperature and monthly temperature and precipitation for 1981–2020 from PRISM. The PRISM web service allows each file to be downloaded twice per day from the same IP address, and the script skips files already on disk, so it can be re-run after an interruption. The second computes the monthly temperature exposure bins (about a day of computing, 10 GB of disk) and aggregates all variables to counties.

The rebuild reproduces the archived files. The county aggregation matrix is identical to the 2021 one, and the exposure bins recomputed from current PRISM files for January 1981 and August 2020 equal the 2021 values exactly.

To also recompute the land cover shares, delete `data2/PRISM/nlcd_prism_weights.tif` and place an NLCD land cover raster (`.img` or `.tif`) in `data/NLCD`. The file used in 2021, `NLCD_2016_Land_Cover_L48_20190424.img`, is part of the CISER archive. With that file, the recomputation takes about half an hour and reproduces the archived shares: 375 of 7.7 million values differ, by at most 0.006 percentage points (about one 30 m pixel in a 4 km grid cell).

## Re-running one figure

The scripts are divided into blocks marked by comment lines such as

```r
# ==== BLOCK: fig8 DEPS: step ====
...
# ==== END BLOCK ====
```

`make block name=fig8` runs the `setup` block of `2_nonlinear_effects.R`, the step-function regressions the figure depends on, and the code for Figure 8, then prints the path of the saved figure. `make list` lists the blocks. The markers are comments, so running a script from top to bottom is unaffected.

| Figure | File | Block |
|---|---|---|
| 1 | `fig1_map_stations.png` | `fig1` |
| 2 | `fig2_map_grid.png` | `fig2` |
| 3 | `fig3_map_us_stations.png` | `fig3` |
| 4 | `fig4_interpol_prism.png` | `fig4` |
| 5 | `fig5_map_landcover_prism.png` | `fig5` |
| 6 | `fig6_temperature_bins.png` | `fig6` |
| 7 | `fig7_map_above30C.png` | `fig7` |
| 8 | `fig8_step.png` | `fig8` |
| 9 | `fig9_spline.png` | `fig9` |
| 10 | `fig10_poly.png` | `fig10` |
| 11 | `fig11_timevarying.png` | `fig11` |
| 12 | `fig12_timevarying.png` | `fig12` |
| 13 | `fig13_spdep.png` | `fig13` |
| 14 | `fig14_specchart.png` | `fig14` |

## How the update was checked

The 2021 versions of scripts 2–5, run on the archived data with their original packages wherever these still install (lfe, sp, raster), reproduce the published Figures 8–14. The updated scripts match these runs. Point estimates agree to about 1e-11 in relative terms and standard errors to about 1e-10 (1e-7 for the spatial error model). The 72 models of Figure 14 have the same estimates and ordering, and confidence intervals that agree to six significant digits.

Pixel by pixel, the figures produced by the updated code compare with the published ones as follows.

- Figures 2 and 5–14 are identical up to anti-aliasing, except for the legends of Figures 9 and 10 (item 6 below).
- Figures 1, 3 and 4 differ in a few pixels along coastlines and county borders. The Natural Earth coastlines differ slightly from those of rworldmap, and sf draws polygon outlines slightly differently than sp. The data shown are the same.

## Notes on the 2021 code

Updating the code turned up a few places where the code or the published figures depart from the text of the chapter. The updated code reproduces the published figures, so it keeps choices 1–3 and flags them in comments.

1. **Sample of counties in Figures 8–12 and 14.** Scripts 2, 3 and 5 project the county map to an Albers projection and then keep counties whose centroid has `x > -100`. Because `x` is in meters after the projection, this keeps counties east of about 96°W (2,127 counties) where the text says east of the 100th meridian (2,510 counties). Script 4 (Figure 13) does not project the map and uses the 100th meridian.
2. **Growing seasons in Figure 14.** The specification labeled March–August uses April–August (months 4 to 8), and the one labeled April–September uses March–September (months 3 to 9). The baseline model of Figure 14 is therefore estimated over April–August.
3. **Sample in Figures 11 and 12.** The time-varying model uses counties in Illinois, Indiana, Iowa and Ohio. The captions describe the sample as counties east of the 100th meridian.
4. **Captions.** The captions of Figures 9, 10 and 12 say that standard errors are clustered by state. The code clusters them by state and year. The caption of Figure 4 describes panels that differ from those in the figure.
5. **Figure 7.** The archived script used 50-hour color intervals up to 900 hours, while the published figure uses 24-hour intervals up to 744 hours and a taller legend. The updated script reproduces the published figure.
6. **Figures 9 and 10.** The legend of the middle panels, which lists 7 basis columns, was drawn with `ncol = 3.5`, and the published version omits column 7. The updated script shows all seven.
7. **Spatial error model (Figure 13).** The 2021 script selected columns by position after `BMisc::makeBalancedPanel()`. Current versions of BMisc return the columns in a different order, which pairs counties with the wrong neighbors and lowers the spatial error coefficient from 0.80 to 0.43. The updated script selects columns by name.
8. **Objects defined interactively.** Script 1 used an object `brks` (Figure 3) and a file list (Figure 6) that were only defined in the interactive session, and scripts 2 and 5 used data.table without loading it, so these scripts failed in a new R session. These are fixed.

## License

The code is released under the MIT License (see `LICENSE`). The data files attached to the release remain subject to the terms of their providers (see [Data](#data)).

## Contact

Ariel Ortiz-Bobea, Cornell University (ao332@cornell.edu)
