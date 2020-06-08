%Xu li so lieu cho do thi ap suat
data = readtable('P at D diagram.csv'); s_D = data{:,1}; P_D = data{:,2};



% S=S-min(S);
% S=S*39.0362/max(S);
% P=P-min(P);
% P=P*15*10^5/max(P);


% plot(S,P,'k')

%Algorithm to reverse array
% Pdao = zeros(326,1);
% Sdao = zeros(326,1);
% for count=1:163
%     Pdao(count) = P(326+1-count);
%     Sdao(count) = S(326+1-count);
%     Pdao(326+1-count) = P(count);
%     Sdao(326+1-count) = S(count);
% end


%for P_B
% s=s-min(s);
% s=s*39/max(s);
% P=P-min(P);
% P=P*23000/max(P);