version 17
clear all
set more off

args outfile
set obs 400
generate double x1 = (mod(_n, 17) - 8) / 5
generate double x2 = mod(_n, 2)
generate double w = 1 + mod(_n, 7) / 10
generate double y = 1 + 0.6*x1 - 0.25*x2 + 0.15*sin(_n/3) + ///
    (0.35 + 0.08*x2)*cos(_n/11)

quietly qrprocess y x1 x2 [pweight=w], ///
    method(onestep, first(qreg)) ///
    qlow(0.1) qhigh(0.9) qstep(0.01) vce(novar) noprint

matrix B = e(coefmat)
clear
svmat double B, names(col)
generate str12 term = ""
replace term = "x1" in 1
replace term = "x2" in 2
replace term = "_cons" in 3
order term
export delimited using "`outfile'", replace
