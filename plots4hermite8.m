function plots4hermite8()
close all
ndata = 19;
data = zeros(9,10,19);
k = 0;
method = ["dijkstra_fs_ibox_slo","dial_fs_iball_slo","dijkstra_freea1_ibox_slo",...
    "dijkstra_es_ibox_slo","dial_es_iball_slo","FMMgrad_ibox_slo","olim8mp0_ibox_slo"];
metname = ["JMM1-Dijkstra","JMM1-Dial","JMM2-Dijkstra","JMM3-Dijkstra","JMM3-Dial","FMM","OLIM8mp0"];
slo = ['1','p','v','m','g'];
sloname = ["Slowness = 1","Linear speed 1","Linear speed 2","Sinefun","Sloth"];
stmax = [4,5,5,5];
len = zeros(19,1);
len(1:14) = 10;
len(15:19) = 6;
mark = ['<','*','s','d','o','>','x'];
col = [ 0,0.6,0;
        0,0.4,0.2;
        0.77,0,0.07;
        0,0,1;
        0.5,0,1;
        0,0,0;
        0.8,0,0.8];
xlabelname = ["N","N","ErrMax/Umax","ERMS/URMS","N","N"];
ylabelname = ["ErrMax/Umax","ERMS/URMS","CPU time","CPU time","Gradient, ErrMax", "Gradient, ERMS"];
figname = ["Figures/ErrMax_vs_N_slo_","Figures/ERMS_vs_N_slo_","Figures/CPUtime_vs_ErrMax_slo_",...
        "Figures/CPUtime_vs_ERMS_slo_","Figures/GradErrMax_vs_N_slo_","Figures/GradERMS_vs_N_slo_"];
ff = 6;
Nmet = 7;
Nslo = 5;
Ntable = 4;
p = zeros(Nslo,Nmet,Ntable);
I = 1:7;
for islo = 1 : Nslo
    for met = 1 : Nmet
        if ~(islo == 4 && (met == 2 || met == 5))               %islo ~= 4 || met ~=4
            fname = sprintf(strcat('Data/',method(met),slo(islo),'.txt'));
            dd = load(fname);
            %% plot Errmax/Umax vs N
            figure((islo-1)*ff + 1);  hold on; 
            plot(dd(:,1),dd(:,4),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            p(islo,met,1) = lsfit(dd(:,1),dd(:,4),fname,'ErrMax/Umax, N:');
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
            %% plot ERMS/URMS vs N
            figure((islo-1)*ff + 2); hold on; 
            plot(dd(:,1),dd(:,5),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            p(islo,met,2) = lsfit(dd(:,1),dd(:,5),fname,'ERMS/URMS, N:');
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
            %% plot CPU time vs Errmax/Umax
            figure((islo-1)*ff + 3); hold on; 
            plot(dd(:,4),dd(:,8),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
            %% plot CPU time vs ERMS/URMS
            figure((islo-1)*ff + 4); hold on; 
            plot(dd(:,5),dd(:,8),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
            %% plot grad ErrMax vs N
            figure((islo-1)*ff + 5);  hold on; 
            plot(dd(:,1),dd(:,6),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            p(islo,met,3) = lsfit(dd(I,1),dd(I,6),fname,'grad ErrMax, N:');
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
            %% plot grad erms vs N
            figure((islo-1)*ff + 6); hold on; 
            plot(dd(I,1),dd(I,7),'Linewidth',2,'Marker',mark(met),'Markersize',10,'color',col(met,:));
            p(islo,met,4) = lsfit(dd(I,1),dd(I,7),fname,'grad ERMS, N:');
            set(gca,'Xscale','log','YScale','log','Fontsize',20);
        end
    end
    for ifig = 1 : ff
        figure((islo-1)*ff + ifig);
        grid;
        title(sloname(islo),'Fontsize',16);
        xlabel(xlabelname(ifig));
        ylabel(ylabelname(ifig));
        if islo ~= 4
            legend(metname,'Location','southwest');
        else % islo == 4           
            legend(metname([1,3,4,6,7]),'Location','southwest');
        end 
        set(gca,'Xscale','log','YScale','log','Fontsize',16);
        saveas(gcf,strcat(figname(ifig),slo(islo)),'epsc');
    end
end

%% Make a preparation for tables
for itable = 1 : Ntable    
    fprintf('\n\n');
    fprintf(strcat('Method &',metname(1),'&',metname(2),'&',metname(3),'&',metname(4),'&',metname(5),'&',metname(6),'&',metname(7),'\\'));
    fprintf('\n');
    for islo = 1 : 5
        fprintf(strcat(sloname(islo),'&'));
        fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\ \n',p(islo,1,itable),...
            p(islo,2,itable),p(islo,3,itable),p(islo,4,itable),p(islo,5,itable),p(islo,6,itable),p(islo,7,itable));
    end
end
end
%%
function rate = lsfit(dx,dy,str1,str2)
    x = log(dx);
    y = log(dy);
    p = polyfit(x,y,1);
    rate = -p(1);
    fprintf(strcat(str1,', ',str2,sprintf(' p = %d, C = %d',rate, exp(p(2)))));
    fprintf('\n');
end
