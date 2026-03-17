Adding a fecundity-survival trade-off to a discrete population model with maturation delay
========

#### Authors
Christopher J. Greyson-Gaito<sup>*1</sup>, Sabrina H. Streipert<sup>2</sup>, Gail S.K. Wolkowicz<sup>1</sup>
---------

### Affiliations
*Corresponding Author - christopher@greyson-gaito.com

1. Department of Mathematics and Statistics, McMaster University, Hamilton, ON, Canada
2. Department of Mathematics, University of Pittsburgh, Pittsburgh, PA, USA

## ORCID
* CJGG &ndash; 0000-0001-8716-0290
* SHS &ndash; 0000-0002-5380-8818
* GSKW &ndash; 0000-0002-4501-2342

[![DOI](https://zenodo.org/badge/916301957.svg)](https://doi.org/10.5281/zenodo.19068255)

## Julia and XPPAUT scripts and datasets

### Folder and file structure
* data &ndash; empty folder for data files to be placed (created in TradeOffs_figurecreation.jl)
* figs &ndash; empty folder for figures to be placed (created in TradeOffs_figurecreation.jl)
* src
    * julia
        * packages.jl &ndash; list of packages required (file used in other scripts)
        * TradeOffs_CommonCode.jl &ndash; julia script file containing common code used in other scripts
        * TradeOffs_figurecreation.jl &ndash; julia script to produce the figures in the manuscript
    * xppaut
        * RickerConstanttau1_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=1)
        * RickerConstanttau1.dat &ndash; data output from RickerConstanttau1_bifdiagram.ode
        * RickerConstanttau2_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=2)
        * RickerConstanttau2.dat &ndash; data output from RickerConstanttau2_bifdiagram.ode
        * RickerConstanttau3_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=3)
        * RickerConstanttau3.dat &ndash; data output from RickerConstanttau3_bifdiagram.ode
        * RickerConstanttau4_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=4)
        * RickerConstanttau4.dat &ndash; data output from RickerConstanttau4_bifdiagram.ode
        * RickerConstanttau5_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=5)
        * RickerConstanttau5.dat &ndash; data output from RickerConstanttau5_bifdiagram.ode
        * RickerConstanttau6_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=6)
        * RickerConstanttau6.dat &ndash; data output from RickerConstanttau6_bifdiagram.ode
        * RickerConstanttau8_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for RickerConstant model (with tau=8)
        * RickerConstanttau8.dat &ndash; data output from RickerConstanttau8_bifdiagram.ode
        * RickerRickertau1_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=1)
        * RickerRickertau1_C.dat &ndash; data output from RickerRickertau1_bifdiagram.ode
        * RickerRickertau1_beta_b.dat &ndash; data output from RickerRickertau1_bifdiagram.ode
        * RickerRickertau2_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=2)
        * RickerRickertau2_C.dat &ndash; data output from RickerRickertau2_bifdiagram.ode
        * RickerRickertau2_beta_b.dat &ndash; data output from RickerRickertau2_bifdiagram.ode
        * RickerRickertau3_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=3)
        * RickerRickertau3_C.dat &ndash; data output from RickerRickertau3_bifdiagram.ode
        * RickerRickertau3_beta_b.dat &ndash; data output from RickerRickertau3_bifdiagram.ode
        * RickerRickertau4_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=4)
        * RickerRickertau4_C.dat &ndash; data output from RickerRickertau4_bifdiagram.ode
        * RickerRickertau4_beta_b.dat &ndash; data output from RickerRickertau4_bifdiagram.ode
        * RickerRickertau5_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=5)
        * RickerRickertau5_C.dat &ndash; data output from RickerRickertau5_bifdiagram.ode
        * RickerRickertau5_beta_b.dat &ndash; data output from RickerRickertau5_bifdiagram.ode
        * RickerRickertau9_bifdiagram.ode &ndash; xppaut script to produce bifurcation diagram for the RickerRicker model (with tau=9)
        * RickerRickertau9_C.dat &ndash; data output from RickerRickertau9_bifdiagram.ode
        * RickerRickertau9_beta_b.dat &ndash; data output from RickerRickertau9_bifdiagram.ode
* .gitignore &ndash; file containing files and folders that git should ignore
* LICENSE.txt &ndash; CC by 4.0 License for this repository
* README.md &ndash; this file
* Manifest.toml &ndash; Record of state of packages in the julia environment
* Project.toml &ndash; Gives julia package dependencies
