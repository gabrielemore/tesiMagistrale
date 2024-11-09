%% IMPORTAZIONE DATI

data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_off\Sim_gabriele_pasti_Circadian_off.mat");

%% ESTRAZIONE DATI

%time full simulation
time = data_circadian_off.SimSettings.simDurationMinutes; %[min]

%G(t) - Blood glucose concetration
G=data_circadian_off.Subjects(1).Signals.subjE_Gp_conc  ; %[mg/dl] Plasma glucose concentration

%Q_i(t) - insulin delivery rate in plasma
IOB_basal = data_circadian_off.Subjects(1).Signals.insulin_Normal_Basal_Use_IOB;  %[units]
IOB_bolus = data_circadian_off.Subjects(1).Signals.insulin_Normal_Bolus_Use_IOB;  %[units]
Q_i = IOB_basal + IOB_bolus;%[units]

%Q_g(t) - rate of carbohydrate absorption from the gut
Q_g=data_circadian_off.Subjects(1).Signals.subjE_G_app_rate; %[mg/kg/ min] Rate of appearance of glucose 

%Q_isub(t) - rate of carbohydrate absorption from the gut
Q_isub=data_circadian_off.Subjects(1).Signals.subjE_Gs_conc; %[mg/dl] Subcutaneous glucose concentration 

%Q_sto(t) - rate of carbohydrate absorption from the gut
Q_sto=data_circadian_off.Subjects(1).Signals.subjE_G_util_rate; %[mg/kg/ min] Glucose utilization rate at all tissues  

%u(t) - insulin infusion
inputs_subqInsulin_Normal_Basal_Use= data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Basal_Use; %[pmol] amount of Normal Basal Insulin injected subcutaneously within the current iteration of the simulation   
inputs_subqInsulin_Normal_Bolus_Use=data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Bolus_Use;  %[pmol] amount of Normal Bolus Insulin injected subcutaneously within the current iteration of the simulation  
u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use; %[pmol] 

%r(t) - rate of carbohydrate(CHO) intake
r=data_circadian_off.Subjects(1).Signals.inputs_mealCHO; %[pmol] meal carbohydrates within the current iteration of the simulation 


%% CONVERSIONE UNITA DI MISURA

%%