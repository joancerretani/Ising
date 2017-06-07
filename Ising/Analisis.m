% Correlation:

Ec = load('Corelation_E.txt');
Mc = load('Corelation_M.txt');

ax1 = subplot(2,1,1);
plot(ax1,Ec(:,1),Ec(:,2))
title('Correlacion E')
xlabel('tau')
ylabel('R(tau)')
xlim([0 max(Ec(:,1))])

ax2 = subplot(2,1,2);
plot(ax2,Mc(:,1),Mc(:,2))
title('Correlacion M')
xlabel('tau')
ylabel('R(tau)')
xlim([0 max(Mc(:,1))])

% Termalization:

Et = load('E_termalizacion.txt');
Mt = load('M_termalizacion.txt');

ax1 = subplot(2,1,1);
plot(ax1,Et(:,1),Et(:,2))
title('Termalizacion E')
xlabel('t')
ylabel('E')

ax2 = subplot(2,1,2);
plot(ax2,Mt(:,1),Mt(:,2))
title('Termalizacion M')
xlabel('t')
ylabel('M')

% M vs T:

MvsT = load('MvsT.txt');
plot(MvsT(:,1),MvsT(:,2))
title('M vs T')
xlabel('T')
ylabel('M')


% Hist M:

Hist_M = load('Hist_M.txt');
time_M = load('time_M.txt');

ax1 = subplot(2,1,1);
plot(ax1,Hist_M(:,1),Hist_M(:,2))
title('hist M')
xlabel('M')
ylabel('Hist')

ax2 = subplot(2,1,2);
plot(ax2,time_M(:,1),time_M(:,2))
title('tiempo M')
xlabel('t')
ylabel('M')

% E cor vs T
corE_T = load('Cor_E_vs_T.txt');
plot(corE_T(:,1),corE_T(:,2),'.')

% Cv
varEvsT = load('varMvsT.txt');
plot(varEvsT(:,1),sqrt(varEvsT(:,2)),'.')
