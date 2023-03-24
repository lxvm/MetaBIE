#!/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
METABIEPATH=$SCRIPTPATH/MetaBIE

julia --project=. -e "import Pkg; Pkg.develop(path=\"$METABIEPATH\"); Pkg.instantiate()"
julia --project=. --threads 1 code/make_plots.jl

mv -t report fields.png scaling.png densities.png params.png densities_conv_direct.png fields_conv_direct.png densities_conv_fast.png fields_conv_fast.png

cd report
latexmk report_VanMunoz_Lorenzo.tex
cd -