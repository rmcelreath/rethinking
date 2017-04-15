
#' Human forager hunting returns data
#'
#' Hunting returns of individual Ache men, 1981 to 2007.
#'
#' @docType data
#' @name Achehunting
#' @usage data(Achehunting)
#' @format \enumerate{
#'   \item month : Month of record
#'   \item day : Day of record
#'   \item year : Year of record
#'   \item id : Identifier of individual man
#'   \item age : Man's age at time of record
#'   \item kg.meat : Kilograms of meat returned from hunt
#'   \item hours : Duration in hours of hunting trip
#'   \item datatype : 1
#' if duration of trip known, 3 otherwise }
#' @source Hill and Kintigh. 2009. Current Anthropology 50:369-377.
NULL





#' Bangladesh contraceptive use data
#'
#' Contraceptive use data from 1934 Bangladeshi women.
#'
#' @docType data
#' @name bangladesh
#' @usage data(bangladesh)
#' @format \enumerate{
#'   \item woman : ID number for each woman in sample
#'   \item district : Number for each district
#'   \item use.contraception : 0/1 indicator of contraceptive use
#'   \item living.children : Number of living children
#'   \item age.centered : Centered age
#'   \item urban : 0/1 indicator of urban context }
#' @source Bangladesh Fertility Survey, 1989
NULL





#' Chimpanzee prosocialty experiment data
#'
#' Data from behavior trials in a captive group of chimpanzees, housed in
#' Lousiana. From Silk et al. 2005. Nature 437:1357-1359.
#'
#' @docType data
#' @name chimpanzees
#' @usage data(chimpanzees)
#' @format \enumerate{
#'   \item actor : name of actor
#'   \item recipient : name of recipient (NA for partner absent condition)
#'   \item condition : partner absent (0), partner present (1)
#'   \item block : block of trials (each actor x each recipient 1 time)
#'   \item trial : trial number (by chimp = ordinal sequence of
#'     trials for each chimp, ranges from 1-72; partner present trials were
#'     interspersed with partner absent trials)
#'   \item prosocial_left : 1 if prosocial (1/1) option was on left
#'   \item chose_prosoc : choice chimp made (0 = 1/0 option, 1 = 1/1 option)
#'   \item pulled_left : which side did chimp pull (1 = left, 0 = right) }
#' @author Richard McElreath
#' @source Silk et al. 2005. Nature 437:1357-1359.
#'
NULL





#' Capuchin monkey contests
#'
#' Data on features and outcomes of intergroup contests among territorial
#' groups of capuchin monkey (Cebus capucinus). Each row contains a single
#' contest between two groups.
#'
#'
#' @docType data
#' @name Crofoot
#' @usage data(Crofoot)
#' @format \enumerate{
#'   \item focal : ID of focal group
#'   \item other : ID of other group
#'   \item dyad : ID of specific dyad pair of groups
#'   \item win : 1 if focal won contest, 0 if other won
#'   \item dist_focal : Distance in meters of focal group from
#'     the center of its home range
#'   \item dist_other : Distance in
#' meters of other group from the center of its home range
#'   \item n_focal : Number of individuals in focal group
#'   \item n_other : Number of individuals in other group
#'   \item m_focal : Number of males in focal group
#'   \item m_other : Number of males in other group
#'   \item f_focal : Number of females in focal group
#'   \item f_other : Number of females in other group }
#' @source M.C. Crofoot, I.C. Gilby, M.C. Wikelski, and R.W. Kays. 2008.
#' PNAS 105:577--581.
#'
NULL




#' Dinosaur growth data
#'
#' Mass by age data for 6 species of dinosaur.
#'
#'
#' @docType data
#' @name Dinasaurs
#' @usage data(Dinasaurs)
#' @format \enumerate{
#'   \item age : Estimated age of specimen at death, in years
#'   \item mass : Estimated body mass at death, in kilograms
#'   \item species : Species name
#'   \item sp_id : Identification number of species }
#' @source Erickson, Rogers, Yerby. 2001. Dinosaurian growth patterns and
#' rapid avian growth rates. Nature 412:429-433.
#'
NULL





#' Dissertation page length data
#'
#' Page lengths of 3037 dissertations filed between 2006 and 2014 at the
#' University of Minnesota.
#'
#'
#' @docType data
#' @name Dissertations
#' @usage data(Dissertations)
#' @format \enumerate{
#'   \item pages: number of pages
#'   \item major: department/program of dissertation
#'   \item year: year of filing }
#' @source Data originally parsed by...
#'
NULL








#' Fishing data
#'
#' Fishing data from visitors to a national park.
#'
#' @docType data
#' @name Fish
#' @usage data(Fish)
#' @format \enumerate{
#'   \item fish_caught : Number of fish caught during visit
#'   \item livebait : Whether or not group used livebait to fish
#'   \item camper : Whether or not group had a camper
#'   \item persons : Number of adults in group
#'   \item child : Number of children in group
#'   \item hours : Number of hours group spent in park }
#'
NULL





#' Urban foxes
#'
#' Data on urban fox groups in England.
#'
#'
#' @docType data
#' @name foxes
#' @usage data(foxes)
#' @format \enumerate{
#'   \item group : ID of group
#'   \item avgfood : Average available food in group's territory
#'   \item groupsize : Size of each group
#'   \item area : Area of group territory
#'   \item weight : Body weight of individual fox }
#' @source Modified from an example in Grafen and Hails.
#'
NULL





#' Prairie dog dispersal data
#'
#' Dispersal and kin residence data for three species of prairie dog, from 1976
#' to 2004. Each row is an individual dispersal record, with associated
#' descriptors.
#'
#'
#' @docType data
#' @name Hoogland
#' @usage data(Hoogland)
#' @format \enumerate{
#'   \item Species : 1 = black-tailed, 2 = Gunnison's, 3 = Utah
#'   \item Year : Year of case record
#'   \item Male : 1 = male, 0 = female
#'   \item NoDisperse : 1 = did not disperse within 12 months of weaning, 0 = dispersed
#'   \item Mother : 1 = mother present for 12 months after weaning, 0 = absent
#'   \item Sisters : Minimal number of littermate sisters present at
#'      dispersal/non-dispersal
#'   \item Bros : Minimal number of littermate brothers present
#'   \item ClanSize : Minimal number of adults present in territory at
#'      time of dispersal/non-dispersal. Includes close kin, distant kin,
#'      immigrants, and focal individual
#'   \item AllKin : Minimal number of mother and
#'      littermates in territory at time of dispersal/non-dispersal }
#' @source Hoogland 2013. Science 339:1205--1207.
#'
NULL





#' Hurricane fatalities and gender of names
#'
#' Data used in Jung et al 2014 analysis of effect of gender of name on
#' hurricane fatalities. Note that hurricanes Katrina (2005) and Audrey (1957)
#' were removed from the data.
#'
#'
#' @docType data
#' @name Hurricanes
#' @usage data(Hurricanes)
#' @format \enumerate{
#'   \item name : Given name of hurricane
#'   \item year : Year of hurricane
#'   \item deaths : number of deaths
#'   \item category : Severity code for storm
#'   \item min_pressure : Minimum pressure, a measure of storm strength;
#'      low is stronger
#'   \item damage_norm : Normalized estimate of damage in dollars
#'   \item female : Indicator variable for female name
#'   \item femininity : 1-11 scale from totally masculine (1)
#'      to totally feminine (11) for name. Average of 9 scores from 9 raters.  }
#' @source Jung et al. 2014. Female hurricanes are deadlier than male
#' hurricanes. PNAS.
#'
#'
NULL





#' Oceanic islands data
#'
#' Data on historic tool complexity and demography in various Oceanic islands
#' societies. There are three data tables: (1) \code{Kline} is the basic data
#' table; (2) \code{Kline2} contains latitude and longitude columns; and (3)
#' \code{islandsDistMatrix} is a matrix of pairwise distances between islands,
#' in thousands of kilometers.
#'
#'
#' @aliases Kline Kline2 islandsDistMatrix islands
#' @docType data
#' @name Kline
#' @usage data(Kline)
#' @format \enumerate{
#'   \item culture : Name of island culture
#'   \item population : Historical population size
#'   \item contact : low or high contact rate with other islands
#'   \item total_tools : number of tools in historical tool kit
#'   \item mean_TU : another measure of tool complexity
#'   \item lat : latitude of island
#'   \item lon : longitude of island
#'   \item lon2 : longitude coded for ease of plotting
#'   \item logpop : natural logarithm of population }
#' @source Kline, M.A. and R. Boyd. 2010. Proc R Soc B 277:2559--2564.
NULL


#' Gift exchange data
#'
#' Data on household gift exchanges from Koster and Leckie. There are two data
#' frames: \code{kl_dyads} and \code{kl_households}.
#'
#' @docType data
#' @name KosterLeckie
#' @usage data(KosterLeckie)
#' @format \code{kl_dyads}: \enumerate{
#'
#'   \item hidA : household ID for A in dyad
#'   \item hidB : household ID for B in dyad
#'   \item did : dyad ID number
#'   \item giftsAB : count of gifts from A to B
#'   \item giftsBA : count of gifts from B to A
#'   \item offset : relative rate of contact in dyad
#'   \item drel1 :
#'   \item drel2 :
#'   \item drel3 :
#'   \item drel4 :
#'   \item dlndist :
#'   \item dass :
#'   \item d0125 :
#'   }
#'   \code{kl_households}: \enumerate{
#'   \item hid : household ID
#'   \item hgame :
#'   \item hfish :
#'   \item hpigs :
#'   \item hwealth :
#'   \item hpastor : }
#' @source Koster and Leckie (2014) Food sharing networks in lowland
#' Nicaragua: An application of the social relations model to count data.
#' Social Networks.
#'
NULL


#' Predictions for map and map2stan models
#'
#' Primate milk data
#'
#' Comparative primate milk composition data, from Table 2 of Hinde and
#' Milligan. 2011. Evolutionary Anthropology 20:9-23.
#'
#'
#' @docType data
#' @name milk
#' @usage data(milk)
#' @format Returns a data frame containing 29 rows and 8 columns.  \enumerate{
#'
#'   \item clade: Broad taxonomic group
#'   \item species: Species name
#'   \item kcal.per.g: Kilocalories per gram of milk
#'   \item perc.fat: Percent fat
#'   \item perc.protein: Percent protein
#'   \item perc.lactose: Percent lactose
#'   \item mass: Body mass of mother, in kilograms
#'   \item neocortex.perc: Percent of brain mass that is neocortex }
#' @author Richard McElreath
#' @source Hinde and Milligan. 2011. Evolutionary Anthropology 20:9-23.
#'
#'
NULL








#' Livestock mortality data
#'
#' Mortality counts for 11 herds of large and small livestock. Large livestock
#' are Zebu cattle, and small livestock are a mix of (fat-tailed, Maasai) sheep
#' and goats.
#'
#'
#' @docType data
#' @name Rinder
#' @usage data(Rinder)
#' @format \enumerate{
#'   \item herd : identifier for individual herd
#'   \item stock : categorical small or large stock
#'   \item n : number of animals at beginning of observation period
#'   \item mortality : number of animals that died in observation period }
#'
NULL



#' Experimental data on ethical dilemmas
#'
#' These data comprise outcomes of experimental ethical dilemmas that are often
#' called 'trolley' problems. Data kindly provided by Fiery Cushman.
#'
#'
#' @docType data
#' @name Trolly
#' @usage data(Trolly)
#' @format \enumerate{
#'   \item case: a code that combines treatment and story labels
#'   \item response: participant's rating of appropriateness of action in
#'      story
#'   \item order: integer order story appeared, within participant
#'   \item id: participant id (factor)
#'   \item age: participant's age in years
#'    \item male: participant's gender; 1 for male, 0 for female
#'   \item edu: participant's highest educational level
#'   \item action: treatment code for
#' story with action (1) or inaction (0)
#'   \item intention: treatment code for intent (1) or lack of intent (0)
#'   \item contact: treatmetn code for contact action (1) or lack of contact (0)
#'   \item story: factor label for basic scenario modified by treatments
#'   \item action2: alternative coding of action that is union of action and
#'      contact variables }
#' @source Cushman et al. 2006. Psychological Science 17:1082--1089.
#'
NULL





#' Ultimate Fighting Championship matches
#'
#' Outcomes of televised Ultimate Fighting Championship (UFC) matches, as a
#' function of handedness of fighters. Data coded by Pollet et al (2013).
#'
#' @docType data
#' @name UFClefties
#' @usage data(UFClefties)
#' @format \enumerate{
#'   \item fight : Unique identifier for match
#'   \item episode : Identifier for UFC episode
#'   \item fight.in.episiode : Order of fight in episode
#'   \item fighter1.win : 1 if fighter 1 won the match; 0 if fight 2 won
#'   \item fighter1 : Unique identifier for fighter 1
#'   \item fighter2 : Unique identifier for fighter 2
#'   \item fighter1.lefty : 1 if fighter 1 was left handed; 0 otherwise
#'   \item fighter2.lefty : 1 if fighter 2 was left handed; 0 otherwise }
#' @source Pollet et al. 2013. Animal Behaviour 86:839--843.
NULL

#' Waffle House and marriage statistics
#'
#' Data for the individual States of the United States, describing number of
#' Waffle House diners and various marriage and demographic facts.
#'
#'
#' @docType data
#' @name WaffleDivorce
#' @usage data(WaffleDivorce)
#' @format \enumerate{
#'   \item Location : State name
#'   \item Loc : State abbreviation
#'   \item Population : 2010 population in millions
#'   \item MedianAgeMarriage: 2005-2010 median age at marriage
#'   \item Marriage : 2009 marriage rate per 1000 adults
#'   \item Marriage.SE : Standard error of rate
#'   \item Divorce : 2009 divorce rate per 1000 adults
#'   \item Divorce.SE : Standard error of rate
#'   \item WaffleHouses : Number of diners
#'   \item South : 1 indicates Southern State
#'   \item Slaves1860 : Number of slaves in 1860 census
#'   \item Population1860 : Population from 1860 census
#'   \item PropSlaves1860 : Proportion of total population that were slaves in 1860 }
#' @source 1860 census data from http://mapserver.lib.virginia.edu.
#' Marriage and divorce rates from 2009 American Community Survey (ACS). Waffle
#' House density data from wafflehouse.com (retrieved January 2012).
NULL


#' Wine tasting scores, Princeton NJ 2012
#'
#' Numerical scores from 9 judges, both French and American, testing 20
#' different wines, both French and American. Judges were blind to the identity
#' of each wine during the testing.
#'
#' @docType data
#' @name Wines2012
#' @usage data(Wines2012)
#' @format \enumerate{
#'   \item judge : Name of judge
#'   \item flight : white or red group of 10 wines
#'   \item wine : Wine ID
#'   \item score : Numerical score of wine
#'   \item wine.amer : 1 for American wines from New Jersey
#'   \item judge.amer : 1 for American judges }
#' @source Raw data from http://www.liquidasset.com/report161.html
NULL

