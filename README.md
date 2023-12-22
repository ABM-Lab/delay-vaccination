Optimal dosing interval:


Program in matlab associated with paper "Deviation from the recommended schedule: Optimal dosing interval for a two-dose vacciantion programme". 

Model details:
A delay-differential model to describe the dynamics of disease spread with a two-dose vaccination. The model incorporated variables such as waning vaccine-induced and naturally-acquired immunity, as well as the capacity for vaccine distribution. We simulated the model and determined the optimal dosing interval as a function of vaccien efficacy against infection and outcomes by comparing the disease burden in delayed vaccination scenarios to that of following the recommended schedule. 


Repository contents and how to use them:

There are five m files and eight data files. 
RRtau_vacsR: run directly and save a .mat file with total number of incidence/hosp/death of primary infection during the DSD (delayed second dose) period when the second dose is administrated in a recommended schedule;
RRtau_vacsD: run directly and save a .mat file with total number of incidence/hosp/death of primary infection during the DSD period when the second dose is adminstrited in a delayed schedule;
RR100_vacsR: run directly and save a .mat file with total number of incidence/hosp/death of primary infection during the first 100 days when the second dose is administrated in a recommended schedule;
RR100_vacsD: run directly and save a .mat file with total number of incidence/hosp/death of primary infection during the first 100 days when the second dose is adminstrited in a delayed schedule;

plot_RR: load different sets of data, and produce the relative reproduction plot during the DSD or first 100 days;

R1_400: data file for inci/hosp/death during DSD period when second dose is in a recommended schedule with R0 = 1.1; 
D1_400: data file for inci/hosp/death during DSD period when second dose is in a delayed schedule with R0 = 1.1;
RR100_R1400: data file for inci/hosp/death during first 100 days when second dose is in a recommended schedule with R0 = 1.1; 
RR100_D1400: data file for inci/hosp/death during first 100 days when second dose is in a delayed schedule with R0 = 1.1; 

R1_400_R18: data file for inci/hosp/death during DSD period when second dose is in a recommended schedule with R0 = 1.8; 
D1_400_R18: data file for inci/hosp/death during DSD period when second dose is in a delayed schedule with R0 = 1.8; 
RR100_R1400_R18: data file for inci/hosp/death during first 100 days when second dose is in a recommended schedule with R0 = 1.8; 
RR100_D1400_R18: data file for inci/hosp/death during first 100 days when second dose is in a delayed schedule with R0 = 1.8; 

Published Studies:
Deviation from the recommended schedule: Optimal dosing interval for a two-dose vacciantion programme, Zhen Wang, Gergely Rost, and Seyed M. Moghadas, 2023, Royal Society Open Science (in Review).






