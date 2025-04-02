% This "primary MATLAB control element plugin" function, called once on each iteration of the simulation, takes the following parameters:
% subjObject - A 1x1 struct object containing subject parameters
%   name:  Subject name, as displayed in the subject selection tab (e.g. adult#001).
%   type1:  Boolean indicating true if the subject is type 1
%   CR:  Carbohydrate ratio (g/Unit)
%   CF:  Correction Factor (mg/dL/Unit)
%   Gb:  Basal (fasting) glucose concentration (mg/dL)
%   BW:  bodyweight (kg)
%   dailyBasalInsulin:  Daily amount of basal "normal" insulin needed to maintain the basal fasting glucose concentration, Gb (Units)
%   OGTT:  The subject's 2-hour OGTT result -- only applicable to T2 and prediabetic subjects (mg/dL)
%
% sensorSigArray - A Nx1 array containing provided sensors
% The order of signals in this array matches the order as entered via the "Edit Input Signals" dialog box within the configuration of this control.
%
% nextMealObject - A 1x1 struct object containing the properties of next scheduled meal
%   amountMg: The amount of carbohydrates in meal before any multiplier is applied.
%   durationInMinutes: The total duration of meal consumption.
%   minutesUntilNextMeal: The remaining minutes before the next meal.
%   bolusMultiplier: The scaling factor to be applied to the size of the meal.
%
% nextExerciseObject - A 1x1 struct object containing the properties of next scheduled exercise
%   intensityFrac: A value from 0.0 to 1.0:
%       0: No exercise
%       0.25:   Light exercise
%       0.5:    Moderate exercise
%       0.655:  Intense exercise
%   durationInMinutes: Total duration of the exercise session
%   minutesUntilNextSession: The minutes before the next scheduled exercise.
%
% timeObject - A 1x1 struct object containing the properties of simulation time
%   minutesPastSimStart:  number of minutes since the start of the simulation
%   daysSinceJan1:  number of days since Jan 1 (0 for Jan 1, 1 for Jan 2, etc.) -- can be used for time of year, to handle possible seasonal effects
%   daysSinceMonday:  number of days since Monday (0 for Monday, 1 for Tuesday, etc.) -- can be used for handling weekly schedules
%   minutesPastMidnight:  time of day, in minutes since midnight (0=midnight, 360=6AM, 720=noon, etc.)
%
% modelInputsToModObject - A 1x1 struct object containing the properties of model input parameters
%   mealCarbsMgPerMin:  rate of consumption of meal carbohydrates (mg/minute)
%   fullMealCarbMgExpectedAtStart:  value indicated only at start of meal (0 at all other times), representing total carbs in the meal (mg).
%   highFatMealFlag:  boolean (true or false) value to be set at the same time as fullMealCarbMgExpectedAtStart, indicating whether this meal is high in fat
%   highProteinMealFlag:  boolean (true or false) value to be set at the same time as fullMealCarbMgExpectedAtStart, indicating whether this meal is high in protein
%   glucOrDextIvInjMgPerMin:  glucose or dextrose provided by IV (mg/minute)
%   glucagonSqInjMg:  rate of glucagon being provided subcutaneously (mg/minute)
%   exerciseIntensityAsFrac:  intensity of exercise as a fraction of full intensity (1.0 = full intensity, 0.0 = no exercise)
%   sqInsulinNormalBasal:  subcutaneous delivery rate of normal insulin being used for basal insulin (pmol/minute)
%   ivInsulinNormalBasal:  intraveneous delivery rate of normal insulin being used for basal insulin (pmol/minute)
%   sqInsulinNormalBolus:  subcutaneous delivery rate of normal insulin being used as a bolus (pmol/minute) -- usually the whole bolus is given in 1 minute
%   ivInsulinNormalBolus:  intraveneous delivery rate of normal insulin being used as a bolus (pmol/minute) -- usually the whole bolus is given in 1 minute
%   sqInsulinUltraRapidBolus:  intraveneous delivery rate of ultra-rapid insulin being used as a bolus (pmol/minute) -- usually the whole bolus is given in 1 minute
%   slowRelInsulinStandardLongActing:  delivery rate of the standard long-acting insulin (pmol/minute) -- usually the whole dose is given in 1 minute
%   sqCustomInsulin1:       subcutaneous delivery rate of the 1st user-defined insulin (pmol/minute), only to be used when this insulin is not of a "slow release" variety
%   ivCustomInsulin1:       intraveneous delivery rate of the 1st user-defined insulin (pmol/minute), only to be used when this insulin is not of a "slow release" variety
%   slowRelCustomInsulin1:  delivery rate of the 1st user-defined insulin, to be used when this insulin is of a slow release variety.
%   sqCustomInsulin2:  subcutaneous delivery rate of the 2nd user-defined insulin (pmol/minute), only to be used when this insulin is not of a "slow release" variety
%   ivCustomInsulin2:  intraveneous delivery rate of the 2nd user-defined insulin (pmol/minute), only to be used when this insulin is not of a "slow release" variety
%   slowRelCustomInsulin2:  delivery rate of the 2nd user-defined insulin, to be used when this insulin is of a slow release variety.
%   ...etc.
% The order in which custom insulins are associated with the ...CustomInsulin1, 2, etc. properties will correspond to the order in the "insulin definitions" configuration page.
%   drugDoses: A struct object containing the dose of a defined drug. To access the defined drug dose, call "drugDose.(drugTypeName).(drugInstanceName)".
% inNamedSignals:
%   This is an optional parameter that can be used to receive signals generated from other control elements as output signals (in their
%   outSigArray parameter -- see below)
%
% outSigArray:
%   This is an optional parameter that can be populated with values that will be contributed to the recorded state history files, in the same way
%   that sensor signals are.  When such output signals are to be contributed in this way, the user must also implement the numOutSignals and
%   outSignalName functions.
%   The values in outSigArray will be made available to control elements running after this control element, via their inNamedSignals parameter
%
%
% runStopStatusToModObject - A 1x1 struct object containing fields which can be optionally modified by this
%                            MATLAB plugin to force the simulation to stop, and indicate a reason for stopping.
%                            If left unchanged, these will always indicate that the simulation should not be stopped.
%                            Fields are as follows:
%   stopRun:  boolean (true or false) value indicating whether the simulation should stop after returning from the runIteration function
%   error:    boolean value indicating whether an error condition is responsible for stopping the simulation
%   message:  text message (as a MATLAB character array) giving an explanation of why the simulation is to be stopped.
%
%
% getsimsettings(fieldname) - The simulation setting parameters can be obtained by invoking the global "getsimsettings()" function. If the function is called without any input argument, getsimsettings() returns a 1x1 struct object containing available simulation setting information in each field. When using with getsimsettings(fieldname), the function returns the value of the designated field "fieldname". The available fields includes: 
%   simDurationMinutes:  the planned duration of the simulation, in minutes
%   subjectName:  Subject name, as displayed in the subject selection tab (e.g. adult#001).
%   populationName:  name of the population to which this subject belongs, as displayed in the subject selection tab (e.g. Type1 Adult).
%   useUniversalRandomSeed:  boolean, with a value of true if the "Use random seed common to all elements" setting is checked on the "general" tab,
%                            indicating that the same initial seed should be used in random number generation across all sensors, control elements,
%                            and delivery elements.
%   universalRandomSeed:     When "useUniversalRandomSeed" is true, the value of the seed to be used with all elements that use random number
%                            generation.  When useUniversalRandomSeed is true, this value should be used by the MATLAB plugin to initialize
%                            random number generation (typically done on the first iteration of this primary function).
%   simName:                 The name of the simulation, as configured in the "General" tab of the DMMS.R's GUI.
%%%%%%%%%%%
% Warning %
%%%%%%%%%%%
% Please do not run any command on an active Matlab engine when it's visible.
% The simulation result may be affected by any command entered in the console.
% To enable debug log functionality, declare "global debugLog" in the plugin.

% Simple titration approach to target BG
function [modelInputsToModObject, outSignalArray, runStopStatusToModObject] = samplePluginControlElement_matlab(...
          subjObject, sensorSigArray, nextMealObject, nextExerciseObject, timeObject, modelInputsToModObject, inNamedSignals, runStopStatusToModObject)
    global getsimsettings; % Allow to use "getsimsettings()" function for retrieving the simulation settings
    global debugLog; % Enable debug log
  try  
    populationName = getsimsettings('populationName');

    %VARIABILI PERSISTENT
    %var temporali
    persistent Ts Tsim %deltaT
    %parametri sistema
    %persistent o0 o1 o2 o3 o4 o5;
    %matrici sistema
    %persistent A_d B_u_d B_r_d C_d E_d;
    %costanti MPC
    %persistent N umin umax x1min x1max CHO_UB tau IOB_s IOB_d ymin ymax p_ipo p_iper p_ipo_x1 p_iper_x1 p_IOB;
    %osservatore
    persistent Aaug Baug Caug Q_kf R_kf Pkk x_aug u_aug;
    %IOB bound RK vet
    persistent N IOB_vet rk_in rk_in_full;
    %var. CasADi
    persistent MPC

    persistent v_xf v_uf v_x_obs v_xa v_ua v_pkk v_kk v_VN v_u_in v_y v_exit v_delta_iper v_delta_ipo v_delta_iper_x1 v_delta_ipo_x1 v_delta_IOB v_time_sol theta_ott
        
    
    %inizializzazione alla prima iterazione variabili persistent
    if timeObject.minutesPastSimStart == 0
        %var temporali
        Tsim = getsimsettings('simDurationMinutes') +1;
        Ts=5;
        deltaT=1;
        
        %parametri sistema
        [o0,o1,o2,o3,o4,o5] = initializeParameter(subjObject.name);
        %matrici sistema
        [A_d,B_u_d,B_r_d,C_d,E_d] = initializeSysMatrix(o0,o1,o2,o3,o4,o5,Ts,deltaT);
        
        %costanti MPC
        dailyBasalInsulin = subjObject.dailyBasalInsulin;
        Ub = dailyBasalInsulin/(24*60);
        CR=o2/o3;     
        
        %osservatore
        [Aaug,Baug,Caug,Q_kf,R_kf,Pkk] = initializeODO(A_d,B_u_d,B_r_d,C_d,E_d,Ts);
        x0 = [subjObject.Gb Ub Ub 0 0]';
        x_aug = [x0;0;0];
        u_aug = [0;0;1];

        %casaADi
        [MPC,N,IOB_s,IOB_d] = initializeProblemMPC(A_d,B_u_d,B_r_d,C_d,E_d,o4,Ub,CR);
        
        %IOB bound
        IOB_vet = PGcreate_IOB_vector(Tsim,IOB_s,IOB_d);
        
        %RK vet
        [rk_in,rk_in_full] = create_RK_random(Tsim+(N*Ts));
            
        %vettori per salvataggio dati intermedi
        v_xf=[];
        v_uf=[];
        v_x_obs=[];
        v_xa=[];
        v_ua=[];
        v_pkk=[];
        v_kk=[];
        v_VN=[];
        v_u_in=[];
        v_y=[];
        v_exit=[];
        v_delta_iper=[];
        v_delta_ipo=[]; 
        v_delta_iper_x1=[];
        v_delta_ipo_x1=[];
        v_delta_IOB=[];
        v_time_sol =[];
        theta_ott = [o0,o1,o2,o3,o4,o5];

    end

    y = sensorSigArray(1);
    u_in=0;

    if mod(timeObject.minutesPastSimStart,Ts) == 0
        %OSSERVATORE
        %stima x0 con osservatore
        [xk_obs,Pkk,Kk]=PG_ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf);
        %non considero i disturbi per MPC (prendo solo i primi 5 stati)
        xk_sim = xk_obs(1:5);

        %costruzione vettore rk_sim di lunghezza N (predizione 6h)
        rk_sim = create_RK_Sim(rk_in,timeObject.minutesPastSimStart+1,Ts,N);

        %MPC
        %risolvo il problema di minimizzazione al passo k
        [uf_sol,xf_sol,V_sol,d_ipo,d_iper,d_ipo_x1,d_iper_x1,d_IOB,xa_sol,ua_sol]=MPC(xk_sim,rk_sim,IOB_vet((timeObject.minutesPastSimStart+1):Ts:timeObject.minutesPastSimStart+1+N*Ts));
        %verifico exit MPC
        exit=get_stats(MPC);

        %salvataggio dati intermedio
        %converto risultato simbolico in risultato numerico
        xf_sim=full(xf_sol);
        uf_sim=full(uf_sol);
        xa_sim=full(xa_sol);
        ua_sim=full(ua_sol);
        Vn=full(V_sol);
        delta_iper_sim = full(d_iper);
        delta_ipo_sim = full(d_ipo);
        delta_iper_x1_sim = full(d_iper_x1);
        delta_ipo_x1_sim = full(d_ipo_x1);
        delta_IOB_sim = full(d_IOB);

        if exit.success==1
            % prendo solo il primo elemento di uf ottima e lo applico al sistema
            u_in=uf_sim(1);
        else
            %se non trovo una soluzione ottima
            u_in=0;
        end

        %matrice stati x controllati con MPC su N=72 (6 ore)
        v_xf = [v_xf;xf_sim];
        %matrice controllo u
        v_uf = [v_uf;uf_sim];
        %matrice stati artificiali
        v_xa = [v_xa;xa_sim];
        %matrice controllo artificiale
        v_ua = [v_ua;ua_sim];
        %matrice valore del costo
        v_VN = [v_VN; Vn];
        %matrice soluzione MPC trovata
        v_exit=[v_exit,exit.success];
        %salvo matrici osservatore per vedere se convergono
        v_kk = [v_kk;Kk];
        v_pkk = [v_pkk;Pkk];
        %salvo valori delta calcolati in funzione MPC
        v_delta_iper = [v_delta_iper;delta_iper_sim];
        v_delta_ipo = [v_delta_ipo;delta_ipo_sim];
        v_delta_iper_x1 = [v_delta_iper_x1;delta_iper_x1_sim];
        v_delta_ipo_x1 = [v_delta_ipo_x1;delta_ipo_x1_sim];
        v_delta_IOB = [v_delta_IOB;delta_IOB_sim];
        %salvo valori xk osservati dall'osservatore
        v_x_obs = [v_x_obs xk_obs];
        %salvo tempi risoluzione MPC di ogni iterazione
        v_time_sol = [v_time_sol;exit.t_wall_total];

        %aggiorno valori aug per prossima iterazione MPC
        x_aug = xk_obs;
        u_aug = [u_in;rk_in(timeObject.minutesPastSimStart+1);1];
    end

    %salvo y misurata
    v_y = [v_y;y];
    %salvo u ingresso
    v_u_in=[v_u_in;u_in];

    %conversione misure
    %input insulin U->pmol/minute
    correctionBolusDelivered = u_in*6000;
    debugLog = [debugLog 'Correction bolus: ' num2str(correctionBolusDelivered)];

    %input RK al simulatore 
    rescueCarbsDelivered = rk_in(timeObject.minutesPastSimStart+1)*1000; %g->mg
    fullRescueCarbsDelivered = rk_in_full(timeObject.minutesPastSimStart+1)*1000;   

    %input controllo e pasti al simulatore
    modelInputsToModObject.mealCarbsMgPerMin = rescueCarbsDelivered;
    modelInputsToModObject.fullMealCarbMgExpectedAtStart = fullRescueCarbsDelivered;
    modelInputsToModObject.sqInsulinNormalBolus = correctionBolusDelivered;
    
    outSignalArray = {};
    %salvataggio workspace finale
    if Tsim == timeObject.minutesPastSimStart +1
        debugLog = [debugLog 'Save succeded '];
        %cambia on con off
        filename = sprintf('SIM_LIN_CIRC_ON/%s_dati_simulazione_T%d.mat', subjObject.name,Tsim);
        save(filename, 'v_xf', 'v_uf', 'v_x_obs', 'v_xa', 'v_ua', 'v_pkk', 'v_kk', 'v_VN', 'v_u_in', 'v_y', 'v_exit', ...
            'v_delta_iper', 'v_delta_ipo', 'v_delta_iper_x1', 'v_delta_ipo_x1', 'v_delta_IOB', 'v_time_sol','rk_in_full','rk_in','theta_ott');
    end
   
  catch ME
      debugLog = [debugLog getReport(ME)];
  end
end
