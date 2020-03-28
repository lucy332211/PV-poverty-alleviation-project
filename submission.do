clear all
set more off
set matsize 10000


use "submission", clear

global dep	"GDPpc Savepc disincomecity disincomerural"
global indp "treat1"
global c1	"secondgdpr pubexinr snhour   agacre   student mktindex gdppcprovince"
global c2	"secondgdpr pubexinr lnagacre student mktindex lnsnhour lngdppcprovince"

/**** Summary Statistics ****/
estpost tabstat $dep $indp $c1, s(mean sd p25 p50 p75 count) c(s)
eststo stats_county

esttab stats_county using "results/sum_stats.tex", replace 	///
cells("mean(fmt(%9.3fc) label(Mean)) sd(fmt(%9.3fc) label(S.D.)) p25(fmt(%9.3fc) label(Q25)) p50(fmt(%9.3fc) label(Median)) p75(fmt(%9.3fc) label(Q75)) count(fmt(%9.0fc) label(N))") 	///
substitute(\centering \centering\scriptsize) 				///
nomtitle nonumber noobs label


/**  Correlation Matrix **/
estpost correlate disincomerural $indp $c1, matrix
eststo corr_x
esttab corr_x using "results/correlation.tex", replace  		///
star(* 0.1 ** 0.05 *** 0.01) b(%6.2f) 							///
substitute(\centering \centering\scriptsize) 					///
nomtitle nonumber noobs label unstack compress	nonotes not 	

reg disincomerural $indp $c1
estat vif

/**** Figure *****/
preserve
bysort Year: egen treat_yr =sum(treat2)
bysort Year: egen county_yr=count(id)
duplicates drop Year, force
tsset Year
gen new_treat=D.treat_yr
replace new_treat=treat_yr if Year==2013
graph bar new_treat, over(Year) blabel(bar, format(%4.0f)) ytitle(Number of Counties) scheme(s2mono)
graph export "results\guangfu_yr.emf", replace
restore

preserve
bysort Year: egen treat_yr =sum(treat1)
bysort Year: egen county_yr=count(id)
duplicates drop Year, force
tsset Year
gen new_treat=D.treat_yr
replace new_treat=treat_yr if Year==2013
graph bar new_treat, over(Year) blabel(bar, format(%4.0f)) ytitle(Number of Counties) scheme(s2mono)
graph export "results\guangfu_yr2.emf", replace
restore

/****** Baseline Regression ******/
preserve

gen 	duration=Year-start_yr
replace duration=0 if duration==.
gen treat1Xduration=treat1*duration
label var treat1Xduration	"SEPEP*DURATION"
label var duration			"DURATION"

reghdfe lndisincomerural  	treat1 secondgdpr pubexinr lnagacre student, 		absorb(CountyID Year) 	vce(cluster CountyID)
eststo t1_1

reghdfe lndisincomerural  	treat1 $c2, 			 							absorb(CountyID Year) 	vce(cluster CountyID)
eststo t1_2

reghdfe lndisincomerural  	treat1 treat1Xduration duration $c2, 			 	absorb(CountyID Year) 	vce(cluster CountyID)
eststo t1_3

local tb1 "t1_1 t1_2 t1_3"
    esttab `tb1' using "results/baseline.tex", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			keep  (treat1 treat1Xduration duration $c2)    																 ///
			order (treat1 treat1Xduration duration $c2)    																 ///		
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 
restore


/**** Endogenous treatment effect *******/
preserve

foreach v of varlist lndisincomerural $c1 lnsnhour lngdppcprovince{
xtreg `v', fe
predict `v'r, res
replace `v'=`v'r
}

xi:etregress lndisincomerural $c2 i.Year, 	treat(treat1=lnsnhour) 					vce(cluster CountyID)
eststo t1_1

xi:etregress lndisincomerural $c2 i.Year, 	treat(treat1=lnsnhour lngdppcprovince) 	vce(cluster CountyID)
eststo t1_2

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/endotreat.rtf", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust ll, labels("Observations" "Number of Counties" "Log Likelihood") fmt(%9.0fc %9.0fc %9.2fc))  ///					
			keep  (1.treat1 $c2 lngdppcprovince)    															 		 ///
			order (1.treat1 $c2 lngdppcprovince)    															 		 ///		
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 ///
			title(Endogenous Treatment Effect Model)
restore

/**** Parallel Assumption *****/
preserve

gen guangfu=inlist(CountyType,2)
gen tdif=Year-start_yr
gen policym3=inlist(tdif,-3)
gen policym2=inlist(tdif,-2)
gen policym1=inlist(tdif,-1)
gen policy0 =inlist(tdif,0,-1)
gen policy1 =inlist(tdif,1)
gen policy2 =inlist(tdif,2)
gen policy3 =inlist(tdif,3)

foreach v in policym3 policym2 policym1 policy0 policy1 policy2 policy3{
gen guangfuX`v'=guangfu*`v'
}
label var guangfuXpolicym3	"SEPEP*(Year=-3)"
label var guangfuXpolicym2	"SEPEP*(Year=-2)"
label var guangfuXpolicym1	"SEPEP*(Year=-1)"
label var guangfuXpolicy0	"SEPEP*(Year=0)"
label var guangfuXpolicy1	"SEPEP*(Year=+1)"
label var guangfuXpolicy2	"SEPEP*(Year=+2)"
label var guangfuXpolicy3	"SEPEP*(Year=+3)"

reghdfe lndisincomerural guangfuXpolicy* policy*, 		absorb(CountyID Year) vce(cluster CountyID) noconstant
eststo t1_1

reghdfe lndisincomerural guangfuXpolicy* policy* $c2, 	absorb(CountyID Year) vce(cluster CountyID) noconstant
eststo t1_2

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/parallel.tex", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			keep  (guangfuXpolicy* $c2)    																	 	 	 	 ///	
			order (guangfuXpolicy* $c2)    																	 	     	 ///
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 
			
coefplot, ///
   keep(guangfuXpolicym3 guangfuXpolicym2 guangfuXpolicym1 guangfuXpolicy0 guangfuXpolicy1 guangfuXpolicy2 guangfuXpolicy3)  ///
   coeflabels(guangfuXpolicym3 = "-3"   ///
   guangfuXpolicym2 = "-2"              ///   
   guangfuXpolicym1 = "-1"              ///
   guangfuXpolicy0  = "0"               ///
   guangfuXpolicy1  = "+1"              ///
   guangfuXpolicy2  = "+2"              ///   
   guangfuXpolicy3  = "+3")        		///
   vertical                             ///
   yline(0)                             ///
   ytitle("Rural Income Growth Rate %")           			///
   xtitle("Time Passage Relative to Year of SEPAP Policy") 	///
   addplot(line @b @at)                 ///
   ciopts(recast(rcap))                 ///
   rescale(100)							///
   scheme(s2mono)
graph export "results\parallel.emf", replace
graph drop _all
restore


/**** Sub-sample by region and economic condition *****/
preserve
astile rich	=GDPpc, by(Year) nq(2)

gen time=.
gen coef=.
gen ci_upper=.	
gen ci_lower=.

reghdfe lndisincomerural  	treat1 $c2 if Location=="西",  absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1, t(West)

replace time	=1											in 1
replace coef	= _b[treat] 								in 1
replace ci_upper=coef+invttail(e(df_r),0.025)*_se[treat] 	in 1
replace ci_lower=coef-invttail(e(df_r),0.025)*_se[treat] 	in 1


reghdfe lndisincomerural  	treat1 $c2 if Location=="东",  absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2, t(East)

replace time	=2 											in 2
replace coef	= _b[treat] 								in 2
replace ci_upper=coef+invttail(e(df_r),0.025)*_se[treat] 	in 2
replace ci_lower=coef-invttail(e(df_r),0.025)*_se[treat] 	in 2

reghdfe lndisincomerural  	treat1 $c2 if rich==1, 		 absorb(CountyID Year) vce(cluster CountyID)
eststo t1_3, t(Low)

replace time	=3											in 3
replace coef	= _b[treat]									in 3
replace ci_upper=coef+invttail(e(df_r),0.025)*_se[treat] 	in 3
replace ci_lower=coef-invttail(e(df_r),0.025)*_se[treat] 	in 3

reghdfe lndisincomerural  	treat1 $c2 if rich==2, 		 absorb(CountyID Year) vce(cluster CountyID)
eststo t1_4, t(High)

replace time	=4											in 4
replace coef	= _b[treat]									in 4
replace ci_upper=coef+invttail(e(df_r),0.025)*_se[treat] 	in 4
replace ci_lower=coef-invttail(e(df_r),0.025)*_se[treat] 	in 4

local tb1 "t1_1 t1_2 t1_3 t1_4"
    esttab `tb1' using "results/subsample.tex", replace   																  ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																  ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			keep  (treat1 $c2)    																	 	 				  ///	
			order (treat1 $c2)    																	 	 				  ///
			mgroups("Region" "County GDP Per Capita", pattern(1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///	
			coeflabels(_cons "Constant")  eqlabels("") 																	  ///			
			substitute(\centering \centering\scriptsize)  						     									  ///	
			label nonotes compress mtitles nodepvars nogaps 															  
		
tsset time
twoway scatter time coef, ylab(1 "West" 2 "East" 3 "Low GDP Per Capita" 4 "High GDP Per Capita", angle(0)) || ///
rcap ci_upper ci_lower time, xline(0,lpattern(dash)) horizontal    	///           
legend(order(1 "point estimate"            							///
2 "95% conf. int.") pos(6))   ytitle("") scheme(s2mono)
graph export "results/treat-subsample.emf", replace		
graph close _all	
restore

/**** Alternative Definition *****/
preserve

reghdfe lndisincomerural 	 treat3 $c2 if !inlist(CountyType,2),  		absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1

reghdfe lndisincomerural 	 treat2 $c2,  								absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2

reghdfe lndisincomerural 	 treat1 $c2, 								absorb(CountyID Year) vce(cluster CountyID)
eststo t1_3

local tb1 "t1_1 t1_2 t1_3"
    esttab `tb1' using "results/definition.tex", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			indicate("County FE=county" "Year FE=time" "Region-Year FE=region", labels("Y" "N")) 						 ///	
			keep  (treat3 treat2 treat1 $c2)    																	 	 ///	
			order (treat3 treat2 treat1 $c2)    																	 	 ///
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps noomitted													 
restore

/**** Placebo Test *****/
preserve

reghdfe lntaxpc 		treat1 $c2,  	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1, t(Ln(Tax Per Capita))

reghdfe lnGDPpc  		treat1 $c2,  	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2, t(Ln(GDP Per Capita))

reghdfe lnSavepc 		treat1 $c2,  	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_3, t(Ln(Saving Per Capita))

local tb1 "t1_1 t1_2 t1_3"
    esttab `tb1' using "results/placebo.tex", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			indicate("County FE=county" "Year FE=time", labels("Y" "N")) 							 					 ///	
			keep  (treat1 $c2)    																	 	 				 ///	
			order (treat1 $c2)    																	 	 				 ///
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress mtitles nodepvars nogaps 															 
restore

/***** Use raw income as dependent variable *****/
preserve

reghdfe disincomerural  	treat1 secondgdpr pubexinr lnagacre student, 		absorb(CountyID Year) 	vce(cluster CountyID)
eststo t1_1

reghdfe disincomerural  	treat1 $c2, 			 							absorb(CountyID Year) 	vce(cluster CountyID)
eststo t1_2

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/rawincome.tex", replace   																 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			indicate("County FE=county" "Year FE=time", labels("Y" "N")) 					 		 					 ///	
			keep  (treat1 $c2)    																	 	 				 ///
			order (treat1 $c2)    																	 	 				 ///		
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 ///
			title(Raw Income as Dependent Variable)
restore

/***** Propensity Score Matching *****/
preserve
estpost ttest GDPpc student lnagacre lnsnhour, by(treat1) unequal
estimates store balance_before

set seed 10101
gen double x=uniform()
sort x

logit treat1 $c2
predict pscore, pr
twoway (kdensity pscore if treat1==1, bwidth(0.1)) (kdensity pscore if treat1==0, bwidth(0.01) lpattern(dash)), legend( label(1 "treated") label(2 "control" ) ) xtitle("propensity score") title(Propensity Score Before Matching) scheme(s2mono) 

psmatch2 treat1, pscore(pscore) outcome(lnagacre) caliper(0.3) n(100) com logit
graph save "results\psm-pre.gph", replace
keep if _weight>0 & _weight~=.
twoway (kdensity pscore if treat1==1, bwidth(0.1)) (kdensity pscore if treat1==0, bwidth(0.1)  lpattern(dash)), legend( label(1 "treated") label(2 "control" ) ) xtitle("propensity score") title(Propensity Score After Matching)  scheme(s2mono) 
graph save "results\psm-post.gph", replace

graph combine "results\psm-pre.gph" "results\psm-post.gph", graphregion(color(white)) rows(1) 
graph export "results\psm-pre-post.emf",replace
graph close _all	
label var treat2	"GUANGFU"
erase "results\psm-pre.gph"
erase "results\psm-post.gph"

reghdfe lndisincomerural  	treat1 secondgdpr pubexinr lnagacre student,  	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1

reghdfe lndisincomerural  	treat1 $c2,  				 					absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/psm.tex", replace   																 	 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			indicate("County FE=county" "Year FE=time", labels("Y" "N")) 							 					 ///	
			keep  (treat1 $c2)    																	 	 				 ///
			order (treat1 $c2)    																	 	 				 ///		
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 

estpost ttest GDPpc student lnagacre lnsnhour, by(treat1) unequal
estimates store balance_after

global balance2  "GDPpc AA student AA lnagacre AA lnsnhour AA"
cap matrix drop B C
mat B= J(8,5,.)
estimates restore balance_before
matrix A1=e(mu_1)'
matrix B1=e(mu_2)'
matrix C1=e(b)'
matrix D1=e(t)'

estimates restore balance_after
matrix A2=e(mu_1)'
matrix B2=e(mu_2)'
matrix C2=e(b)'
matrix D2=e(t)'

forv i=1/4{
mat S`i' = J(2,5,.) 
	mat S`i'[1,2] = (A1[`i',1],B1[`i',1],C1[`i',1],D1[`i',1]) 
	mat S`i'[2,2] = (A2[`i',1],B2[`i',1],C2[`i',1],D2[`i',1]) 
	mat S`i'[1,1] =999
	mat S`i'[2,1] =998
mat B = B\S`i'
}
mat C = B[9...,1...]
    mat colnames C = Sample Control Treatment Diff T-stats
	mat rownames C= $balance2
	mat list C, format(%6.2f) nodotz noheader

esttab matrix(C, fmt(%9.2f)) using "results/psm.tex", append 	///
nomtitle nonumber noobs label 																											


/**** Parallel Assumption After PSM *****/
gen guangfu=inlist(CountyType,2)
gen tdif=Year-start_yr
gen policym3=inlist(tdif,-3)
gen policym2=inlist(tdif,-2)
gen policym1=inlist(tdif,-1)
gen policy0 =inlist(tdif,0)
gen policy1 =inlist(tdif,1)
gen policy2 =inlist(tdif,2)
gen policy3 =inlist(tdif,3)

foreach v in policym3 policym2 policym1 policy0 policy1 policy2 policy3{
gen guangfuX`v'=guangfu*`v'
}
label var guangfuXpolicym3	"GUANGFU*(Year=-3)"
label var guangfuXpolicym2	"GUANGFU*(Year=-2)"
label var guangfuXpolicym1	"GUANGFU*(Year=-1)"
label var guangfuXpolicy0	"GUANGFU*(Year=0)"
label var guangfuXpolicy1	"GUANGFU*(Year=+1)"
label var guangfuXpolicy2	"GUANGFU*(Year=+2)"
label var guangfuXpolicy3	"GUANGFU*(Year=+3)"

reghdfe lndisincomerural guangfuXpolicy* policy*, 		absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1

reghdfe lndisincomerural guangfuXpolicy* policy* $c2, 	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/parallel-psm.tex", replace   															 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///					
			keep  (guangfuXpolicy* $c2)    																	 	 	 	 ///	
			order (guangfuXpolicy* $c2)    																	 	     	 ///
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress nomtitles nodepvars nogaps 															 
			
coefplot, ///
   keep(guangfuXpolicym3 guangfuXpolicym2 guangfuXpolicym1 guangfuXpolicy0 guangfuXpolicy1 guangfuXpolicy2 guangfuXpolicy3)  ///
   coeflabels(guangfuXpolicym3 = "-3"   ///
   guangfuXpolicym2 = "-2"              ///   
   guangfuXpolicym1 = "-1"              ///
   guangfuXpolicy0  = "0"               ///
   guangfuXpolicy1  = "+1"              ///
   guangfuXpolicy2  = "+2"              ///   
   guangfuXpolicy3  = "+3")        		///
   vertical                             ///
   yline(0)                             ///
   ytitle("Rural Income Growth Rate %") ///
   xtitle("Time Passage Relative to Year of SEPAP Policy") ///
   addplot(line @b @at)                 ///
   ciopts(recast(rcap))                 ///
   rescale(100)                         ///
   scheme(s2mono)
graph export "results\parallel-psm.emf", replace
graph drop _all
restore

/***** Exclude Poor County ****/
preserve

reghdfe lndisincomerural  	treat1 $c2 							 if !inlist(CountyType,1), 		absorb(CountyID Year) vce(cluster CountyID)
eststo t1_1, t(Exclude Guangfu and Poor Counties)

reghdfe lndisincomerural  	treat1 $c2 							 if !inlist(CountyType,1,3), 	absorb(CountyID Year) vce(cluster CountyID)
eststo t1_2, t(Exclude All Poor Counties)

local tb1 "t1_1 t1_2"
    esttab `tb1' using "results/excludepoor.tex", replace   															 ///
	        star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) t(%6.2f)																 ///
            stats(N N_clust r2_within, labels("Observations" "Number of Counties" "Adjusted $ R^2$") fmt(%9.0fc %9.0fc %9.2fc)) ///
			keep  (treat1 $c2)    																	 	 				 ///
			order (treat1 $c2)    																	 	 				 ///		
			coeflabels(_cons "Constant")  eqlabels("") 																	 ///			
			substitute(\centering \centering\scriptsize)  						     									 ///	
			label nonotes compress mtitles nodepvars nogaps 															 
restore
