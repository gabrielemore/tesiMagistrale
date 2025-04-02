function grafici_cvga(root,data,Ts)

npat = 10;
lenghtSim = 1440;
CGM = zeros(lenghtSim, npat); % POPOLARE QUESTA TABELLA CON I VALORI DEL CGM DEI PAZIENTI IN COLONNA


%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

y_mat=[];

for patient=1:10
    if patient ~= 5
    [~,~,~,~,~,~,~,y,~,~] = data_extraction(data_real,patient);

    file_id = sprintf('adult#%03d_dati_simulazione_T10081.mat', patient);
    % file_id = sprintf('adult#%03d_dati_simulazione_T20161.mat', patient);
    file_name = data_patients + file_id;

    load(file_name);

    % Considero solo una settimana
    y=y(1:end/2);
    % v_x_obs=v_x_obs(:,1:floor(end/2)+1);
    % IOB=IOB(1:end/2);
    % v_u_in=v_u_in(1:floor(end/2));
    % rk_in = rk_in(1:floor(end/2)+1);
    % rk_in_full = rk_in_full(1:floor(end/2)+1);
    % v_exit =v_exit(1:floor(end/2)+1);

    % salvo in una matrice
    y_mat= [y_mat y];
    % IOB_mat = [IOB_mat ;IOB'];
    % u_in_mat = [u_in_mat ;v_u_in'];
    % 
    % %%media per IOB vincolo medio
    % o4_mat=[o4_mat theta_ott(5)];
    % CR_mat=[CR_mat CR];
    % Ub_mat=[Ub_mat Ub];
    end
end

CGM=y_mat;


%% CVGA 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2769756/

% Glucose
A.subjE_Gp=CGM(:, 1)'; 
B.subjE_Gp=CGM(:, 2)';
C.subjE_Gp=CGM(:, 3)';
D.subjE_Gp=CGM(:, 4)';
E.subjE_Gp=CGM(:, 5)';
F.subjE_Gp=CGM(:, 6)'; 
G.subjE_Gp=CGM(:, 7)';
H.subjE_Gp=CGM(:, 8)';
I.subjE_Gp=CGM(:, 9)';
%J.subjE_Gp=CGM(:, 10)';

for i = 1:size(CGM, 2)
    Min_gli_tot(i) = prctile(CGM(:,i),2.5);
    Max_gli_tot(i) = prctile(CGM(:,i),97.5);
end
p=polyfit([110 180 300 400],[0 20 40 60],3);
Max_gli_tot = polyval(p, Max_gli_tot);

figure('Position', [100, 100, 700, 700],'Color', 'white')
hold on 
X_pa = [0 20 0]; Y_pa = [0 0 20]; fill(X_pa,Y_pa,[0 1 0], 'HandleVisibility','off'); % A
X_pbpiu = [20 40 0 0]; Y_pbpiu = [0 0 40 20]; fill(X_pbpiu,Y_pbpiu,[7/255 135/255 0/255], 'HandleVisibility','off'); % B 
X_pcpiu = [40 60 0 0]; Y_pcpiu = [0 0 60 40]; fill(X_pcpiu,Y_pcpiu,[1 1 0], 'HandleVisibility','off'); % C
X_pe = [60 60 0]; Y_pe = [0 60 60]; fill(X_pe,Y_pe,[1 0 0], 'HandleVisibility','off'); % D
% %rectangle('Position',[0 1 2 3],'Curvature',[1 1]) % rettangolo con vertice basso sinistra in 0,1 poi a destra di 2 e alto di 3
                                    % % Curvature [1 1] per avere il cerchio
rectangle('Position',[-80 -80 180 180],'Curvature',[1 1],'FaceColor',[1 0 0], 'HandleVisibility','off') %D
hold on
rectangle('Position',[-60 -60 120 120],'Curvature',[1 1],'FaceColor',[1 1 0], 'HandleVisibility','off') %C
hold on
rectangle('Position',[-40 -40 80 80],'Curvature',[1 1],'FaceColor',[7/255 135/255 0/255], 'HandleVisibility','off') %B
hold on
rectangle('Position',[-20 -20 40 40],'Curvature',[1 1],'FaceColor',[0 1 0], 'HandleVisibility','off') %A 
xlim([0 60])
ylim([0 60])
%
text(8,8,'\textbf{A}','HorizontalAlignment','center','Interpreter','Latex','Fontsize',15);
text(30,10,'\textbf{Lower B}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);
text(10,30,'\textbf{Upper B}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);
text(35,35,'\textbf{C}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);
text(50,10,'\textbf{Lower C}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);
text(10,50,'\textbf{Upper C}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);
text(50,50,'\textbf{D}','HorizontalAlignment','center','Interpreter','Latex' ,'Fontsize',15);

set(gca,'YTick',[0 20 40 60])
set(gca,'XTick',[0 20 40 60])
set(gca,'XTickLabel',{'>110' '90' '70' '<50'})
set(gca,'YTickLabel',{'<110' '180' '300' '>400'})

hold on

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; ...
          0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 1 0 0; 0.8 0 0.8; 0 0 1; 0.5 0.5 1];
for i=1:size(Min_gli_tot, 2)
    x1 = min(max(110 - Min_gli_tot(i),0),60);
    y1 = min(max(Max_gli_tot(i),0),60);
   
    PROGR = scatter(x1,y1,500,'filled','MarkerFaceColor',colors(i, :),'MarkerEdgeColor',[0 0 0]);
    hold on
 end

xlabel('Minimum BG','Interpreter','latex'),ylabel('Maximum BG','Interpreter','latex')
grid on
box on
% title('CVGA','Interpreter','latex')

set(findall(gcf,'-property','FontSize'),'FontSize',25);

% salvataggio grafici
% output_folder = 'C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\export_figures\linear_control';
output_folder = 'C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\export_figures\non_linear_control';
% % Esporta in PDF nella cartella desiderata
% filename = fullfile(output_folder, 'Sim_7days_CVGA.pdf');
% export_fig(filename, '-pdf');