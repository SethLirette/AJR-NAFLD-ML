
/*
cd .
*/

import excel "trivariant analysis - NAFLD Score Datasheet anonymized" , case(lower) firstrow clear
drop if missing(age)

*fixing some weird display issues on graphs
replace bmi = 20.99 if bmi==20.9802571
replace bmi = 27.19 if bmi==27.1940091
replace bmi = 29.38 if bmi==29.383953
replace bmi = 38.26 if bmi==38.2567086
replace bmi = 41.48 if bmi==41.4763545

*visual checks
hist nafld
lincheck nafld age
lincheck nafld bmi
lincheck nafld wccm
lincheck nafld sadcm
pwcorr nafld age bmi wccm sadcm

*OLS
reg nafld age bmi wccm sadcm

*MFP
mfp , select(0.05): reg nafld age bmi wccm sadcm // only age and bmi, both stayed linear

*Random Forest
rforest nafld age bmi wccm sadcm , type(reg)
mat list e(importance)

******************************************************************************
*** Bootstrap Block **********************************************************
******************************************************************************

program define bigboot , rclass

    *OLS
	reg nafld age bmi
		  predict yhat  // predicted values
		  
		  gen bias = nafld - yhat

		  gen res_sq  = bias^2 // used for RMSE
		  
		  gen _p30 = (abs(nafld - yhat) / nafld) < 0.3

		summ res_sq
		  return scalar rmse_rg = sqrt(r(mean)) // RMSE
		  
		summ bias , d
		  return scalar bias_rg = r(p50) // median diff bias
		  return scalar prec_rg = r(p75) - r(p25) // precision
		  
		summ _p30
		  return scalar p30_rg = r(mean) // p30
		  		  
		corr nafld yhat
		  return scalar r2_rg = r(rho)^2 // R2
		  
	drop yhat bias res_sq _p30	  
		  
	*Random Forest	  
     rforest nafld age bmi wccm sadcm , type(reg)
		  predict yhat  // predicted values
		  
		  gen bias = nafld - yhat

		  gen res_sq  = bias^2 // used for RMSE
		  
		  gen _p30 = (abs(nafld - yhat) / nafld) < 0.3

		summ res_sq
		  return scalar rmse_rf = sqrt(r(mean)) // RMSE
		  
		summ bias , d
		  return scalar bias_rf = r(p50) // median diff bias
		  return scalar prec_rf = r(p75) - r(p25) // precision
		  
		summ _p30
		  return scalar p30_rf = r(mean) // p30
		  		  
		corr nafld yhat
		  return scalar r2_rf = r(rho)^2 // R2
		  
	drop yhat bias res_sq _p30	  
	  
	*Random Forest with only BMI and Age
     rforest nafld age bmi  , type(reg)
		  predict yhat  // predicted values
		  
		  gen bias = nafld - yhat

		  gen res_sq  = bias^2 // used for RMSE
		  
		  gen _p30 = (abs(nafld - yhat) / nafld) < 0.3

		summ res_sq
		  return scalar rmse_rf2 = sqrt(r(mean)) // RMSE
		  
		summ bias , d
		  return scalar bias_rf2 = r(p50) // median diff bias
		  return scalar prec_rf2 = r(p75) - r(p25) // precision
		  
		summ _p30
		  return scalar p30_rf2 = r(mean) // p30
		  		  
		corr nafld yhat
		  return scalar r2_rf2 = r(rho)^2 // R2
		  
	drop yhat bias res_sq _p30	  
	  
	
end

*Now let's bootstrap our program
bootstrap ///
     rmse_rg=r(rmse_rg) rmse_rf=r(rmse_rf) rmse_rf2=r(rmse_rf2) ///
     bias_rg=r(bias_rg) bias_rf=r(bias_rf) bias_rf2=r(bias_rf2) ///
	 prec_rg=r(prec_rg) prec_rf=r(prec_rf) prec_rf2=r(prec_rf2) ///
	 p30_rg=r(p30_rg)   p30_rf=r(p30_rf)   p30_rf2=r(p30_rf2) ///
	 r2_rg=r(r2_rg)     r2_rf=r(r2_rf)     r2_rf2=r(r2_rf2) ///
	 , reps(500) : bigboot

******************************************************************************
******************************************************************************


******************************************************************************
*** Some plots  **************************************************************
******************************************************************************


*NAFLD Score vs. predicted values
reg     nafld age bmi 
 predict regpred
rforest nafld age bmi wccm sadcm , type(reg)
 predict rfpred
rforest nafld age bmi  , type(reg)
 predict rfpred2
 
twoway ///
  (lfit    nafld nafld if nafld<5 & nafld>-5, lpattern(dash) lcolor(gray)) ///
  (scatter nafld regpred, msymbol(Oh) mcolor(orange%50)) ///
  (scatter nafld rfpred , msymbol(Oh) mcolor(midgreen%50)) ///
  (scatter nafld rfpred2 , msymbol(Oh) mcolor(purple%50)) ///
  , legend(order(1 "Line of Perfect Prediction" 2 "Regression Model" 3 "Random Forest Model" 4 "Alt Random Forest Model") ///
           ring(0) position(5) cols(1)) ///
    xtitle("Predicted NAFLD Score") ytitle("Actual NAFLD Score") 


*Reg Model Bland-Altman
egen avg = rowmean(nafld regpred)
 gen dif = nafld-regpred
summ dif
 global mn = r(mean)
 global ll = r(mean) - 1.96*r(sd)
 global ul = r(mean) + 1.96*r(sd)
 
twoway ///
  (scatter dif avg, msymbol(Oh) mcolor(orange)) ///
  , legend(off) ///
    xtitle("Mean NAFLD Score & Regression Predicted NAFLD Score") ytitle("NAFLD Score - Regression Predicted NAFLD Score") ///
	yscale(range(-5 5)) ///
	ylabel(-5 -2 0 2 5) ///
	xscale(range(-5 5)) ///
	xlabel(-5 -2 0 2 5) ///
	yline($mn $ll $ul, lpattern(solid) lcolor(orange*.6)) ///
	yline(0, lpattern(dash) lcolor(black)) 

drop avg dif


*Random Forest Model Bland-Altman
egen avg = rowmean(nafld rfpred)
 gen dif = nafld-rfpred
summ dif
 global mn = r(mean)
 global ll = r(mean) - 1.96*r(sd)
 global ul = r(mean) + 1.96*r(sd)
 
twoway ///
  (scatter dif avg, msymbol(Oh) mcolor(midgreen)) ///
  , legend(off) ///
    xtitle("Mean NAFLD Score & Random Forest Predicted NAFLD Score") ytitle("NAFLD Score - Random Forest Predicted NAFLD Score") ///
	yscale(range(-5 5)) ///
	ylabel(-5 -2 0 2 5) ///
	xscale(range(-5 5)) ///
	xlabel(-5 -2 0 2 5) ///
	yline($mn $ll $ul, lpattern(solid) lcolor(midgreen*.6)) ///
	yline(0, lpattern(dash) lcolor(black)) 

drop avg dif


*Random Forest2 Model Bland-Altman
egen avg = rowmean(nafld rfpred2)
 gen dif = nafld-rfpred2
summ dif
 global mn = r(mean)
 global ll = r(mean) - 1.96*r(sd)
 global ul = r(mean) + 1.96*r(sd)
 
twoway ///
  (scatter dif avg, msymbol(Oh) mcolor(purple)) ///
  , legend(off) ///
    xtitle("Mean NAFLD Score & Alt Random Forest Predicted NAFLD Score") ytitle("NAFLD Score - Alt Random Forest Predicted NAFLD Score") ///
	yscale(range(-5 5)) ///
	ylabel(-5 -2 0 2 5) ///
	xscale(range(-5 5)) ///
	xlabel(-5 -2 0 2 5) ///
	yline($mn $ll $ul, lpattern(solid) lcolor(purple*.6)) ///
	yline(0, lpattern(dash) lcolor(black)) 

drop avg dif

/*
summ regpred
 di r(sd)/r(mean)
summ rfpred
 di r(sd)/r(mean)
summ rfpred2
 di r(sd)/r(mean)
*/

gen m=(regpred+nafld)/2
gen s2m2=((regpred-nafld)^2/2)/m^2
summ s2m2
di sqrt(r(mean)) // 6.35
drop m s2m2

gen m=(rfpred+nafld)/2
gen s2m2=((rfpred-nafld)^2/2)/m^2
summ s2m2
di sqrt(r(mean)) // 4.51
drop m s2m2

gen m=(rfpred2+nafld)/2
gen s2m2=((rfpred2-nafld)^2/2)/m^2
summ s2m2
di sqrt(r(mean)) // 1.70
drop m s2m2

******************************************************************************
******************************************************************************




