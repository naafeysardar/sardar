drop _all

* THIS FILE GRAPHS THE SHARE OF US GOODS IMPORTS FROM CHINA 1987 TO 2007 (Autor-Dorn-Hanson, Figure 1)

set scheme s2color

use ../dta/figure1_data.dta, clear

twoway (line oil year, lpattern(solid) yaxis(1)) (line cpsman year, lpattern(dash) yaxis(2)) if year>1992 & year<2008, xlab(1993(2)2007) legend(cols(1) lab(1 "Crude Oil Price: WTI $/barrel") lab(2 "Manufacturing employment/Population"))

twoway (line impr year, lpattern(solid) yaxis(1)) (line oil year, lpattern(dash) yaxis(2)) if year>1992 & year<2008, xlab(1993(2)2007) legend(cols(1) lab(1 "China import penetration ratio") lab(2 "Crude Oil Price: WTI $/barrel"))

twoway (line oildd year, lpattern(solid) yaxis(1)) (line oiladhat year, lpattern(dash) yaxis(2)) if year>1992 & year<2008, xlab(1993(2)2007) legend(cols(1) lab(1 "AD Shock (Emerging Markets)") lab(2 "Precautionary Demand Shock"))

twoway (line gasoline year, lpattern(solid) yaxis(1)) (line motor year, lpattern(dash) yaxis(2)) if year>1992 & year<2008, xlab(1993(2)2007) legend(cols(1) lab(1 "RPCE Gasoline") lab(2 "RPCE Motor Vehicles"))

twoway (line gasoline year, lpattern(solid) yaxis(1)) (line clothing year, lpattern(dash) yaxis(2)) if year>1992 & year<2008, xlab(1993(2)2007) legend(cols(1) lab(1 "RPCE Gasoline") lab(2 "RPCE Clothing"))
