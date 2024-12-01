function [time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian,patient)

%time full simulation [min]
time = data_circadian.SimSettings.simDurationMinutes;

%Gb - Glucose basal rate [mg/dl]
Gb= data_circadian.Subjects(patient).Params.Gb;

%CR - Ratio carbohydrate-to-insuline ratio [g/U]
CR=data_circadian.Subjects(patient).Params.CR;

%CF - Correction factor [mg/(dL·U)]
CF=data_circadian.Subjects(patient).Params.CF;

%Ub_day - daily basal insulin [U/day]
Ub_day = data_circadian.Subjects(patient).Params.dailyBasalInsulin;
% Ub - [U/day]->[U/min]
Ub = Ub_day/(24*60);

% u(t) - Input insuline [pmol]
inputs_subqInsulin_Normal_Basal_Use= data_circadian.Subjects(patient).Signals.inputs_subqInsulin_Normal_Basal_Use(2:end);
inputs_subqInsulin_Normal_Bolus_Use=data_circadian.Subjects(patient).Signals.inputs_subqInsulin_Normal_Bolus_Use(2:end);
u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use;
%[pmol]->[U]
u = u/6000;

%r(t) - rate of carbohydrate(CHO) intake [mg]
r=data_circadian.Subjects(patient).Signals.inputs_mealCHO(2:end);
%[mg]->[g]
r = r/1000;

%y(t) - sensor [mg/dl]
y = data_circadian.Subjects(patient).Signals.Sensors__replace_me__(2:end);

%IOB [U]
IOB_b = data_circadian.Subjects(patient).Signals.insulin_Normal_Basal_Use_IOB(2:end);
IOB_bol = data_circadian.Subjects(patient).Signals.insulin_Normal_Bolus_Use_IOB(2:end);
IOB = IOB_b + IOB_bol;

%Rate of glucose apperance in plasma [mg/kg/ min]
Ra = data_circadian.Subjects(patient).Signals.subjE_G_app_rate(2:end);
Vg = 2.1143; % [dL/kg];
Ra=Ra/Vg;

end

