# Run the R scripts in handbook_R/ from the root of the repository.
#
#   make all              download the data and run scripts 0-5 in order
#   make data             0_download_data.R
#   make weather          1_weather_data.R           (figures 1-7)
#   make nonlinear        2_nonlinear_effects.R      (figures 8-10)
#   make timevarying      3_time-varying_effects.R   (figures 11-12)
#   make spatial          4_spatial_dependence.R     (figure 13)
#   make robustness       5_robustness_checks.R      (figure 14)
#   make block name=fig8  re-run one marked block (e.g. a single figure)
#   make list             list the marked blocks

R := Rscript
CODE := handbook_R

.PHONY: all data weather nonlinear timevarying spatial robustness block list

all: data weather nonlinear timevarying spatial robustness

data:
	cd $(CODE) && $(R) 0_download_data.R

weather:
	cd $(CODE) && $(R) 1_weather_data.R

nonlinear:
	cd $(CODE) && $(R) 2_nonlinear_effects.R

timevarying:
	cd $(CODE) && $(R) 3_time-varying_effects.R

spatial:
	cd $(CODE) && $(R) 4_spatial_dependence.R

robustness:
	cd $(CODE) && $(R) 5_robustness_checks.R

block:
	cd $(CODE) && $(R) _run_block.R $(name)

list:
	@grep -ho '^# ==== BLOCK: [^ ]*' $(CODE)/[0-9]*.R | sed 's/# ==== BLOCK: //' | grep -v '^setup$$'
