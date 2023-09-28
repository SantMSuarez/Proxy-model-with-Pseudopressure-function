
WELLS_PARAM = readtable('Wells_GOR.xlsx')
WELLS_NAMES=WELLS_PARAM.Well_name;
WELLS_qo_Start_date=WELLS_PARAM.qo_Start_date;	
WELLS_qg_Start_date=WELLS_PARAM.qg_Start_date;	
WELLS_qw_Start_date=WELLS_PARAM.qw_Start_date;

WELLS_GOR=num2str(round(WELLS_PARAM.GOR_ft3bbl));

WELLS_R2=zeros(length(WELLS_NAMES),1);
WELLS_RMSE=zeros(length(WELLS_NAMES),1);
WELLS_R2_Ori=zeros(length(WELLS_NAMES),1);
WELLS_RMSE_Ori=zeros(length(WELLS_NAMES),1);

for i=1:length(WELLS_NAMES)

temp_str=strcat('INPUT_PROXY_',WELLS_NAMES(i),'.xlsx');
temp_PROXY=readtable(char(temp_str));
temp_X(:,1)=log(pseudop(temp_PROXY.WHP_psig));
temp_X(:,2)=log(pseudop(temp_PROXY.Pr_psia));
temp_X(:,3)=temp_PROXY.WGR_STB_MMSCF;
temp_X(:,4)=temp_X(:,1).*temp_X(:,2);
temp_X(:,5)=temp_X(:,1).*temp_X(:,3);
temp_X(:,6)=temp_X(:,2).*temp_X(:,3);
temp_X(:,7)=temp_X(:,1).^2;
temp_X(:,8)=temp_X(:,2).^2;
temp_X(:,9)=temp_X(:,3).^2;
temp_X(:,10)=temp_X(:,1).*temp_X(:,2).*temp_X(:,3);
temp_X(:,11)=temp_X(:,1).^2.*temp_X(:,2);
temp_X(:,12)=temp_X(:,1).^2.*temp_X(:,3);
temp_X(:,13)=temp_X(:,2).^2.*temp_X(:,1);
temp_X(:,14)=temp_X(:,2).^2.*temp_X(:,3);
temp_X(:,15)=temp_X(:,3).^2.*temp_X(:,1);
temp_X(:,16)=temp_X(:,3).^2.*temp_X(:,2);
temp_X(:,17)=temp_X(:,1).^0.5;
temp_X(:,18)=temp_X(:,2).^0.5;
temp_X(:,19)=temp_X(:,3).^0.5;

temp_Y(:,1)=temp_PROXY.Qg_MMSCF;

% Calculating inverse matrix to estimate proxy parameters
temp_PARAM=inv(temp_X'*temp_X)*temp_X'*temp_Y;

PROXY_PARAM{i}=temp_PARAM;

% If needed, plot predicted versus real gas production 
figure
plot(temp_Y, temp_X*PROXY_PARAM{i},'*r'), hold on

%check R2 and RMSE
WELLS_R2(i)=1-(sum((temp_Y-(temp_X*PROXY_PARAM{i})).^2))...
    /(sum((temp_Y-(mean(temp_Y))).^2));
WELLS_RMSE(i)=sqrt(mean((temp_Y-(temp_X*PROXY_PARAM{i})).^2));

% Clear calculated temporal vectors
clear temp_*

end
% Original
for i=1:length(WELLS_NAMES)

temp_str=strcat('INPUT_PROXY_',WELLS_NAMES(i),'.xlsx');
temp_PROXY=readtable(char(temp_str));
temp_X(:,1)=(temp_PROXY.WHP_psig);
temp_X(:,2)=(temp_PROXY.Pr_psia);
temp_X(:,3)=temp_PROXY.WGR_STB_MMSCF;
temp_X(:,4)=temp_X(:,1).*temp_X(:,2);
temp_X(:,5)=temp_X(:,1).*temp_X(:,3);
temp_X(:,6)=temp_X(:,2).*temp_X(:,3);
temp_X(:,7)=temp_X(:,1).^2;
temp_X(:,8)=temp_X(:,2).^2;
temp_X(:,9)=temp_X(:,3).^2;
temp_X(:,10)=temp_X(:,1).*temp_X(:,2).*temp_X(:,3);
temp_X(:,11)=temp_X(:,1).^2.*temp_X(:,2);
temp_X(:,12)=temp_X(:,1).^2.*temp_X(:,3);
temp_X(:,13)=temp_X(:,2).^2.*temp_X(:,1);
temp_X(:,14)=temp_X(:,2).^2.*temp_X(:,3);
temp_X(:,15)=temp_X(:,3).^2.*temp_X(:,1);
temp_X(:,16)=temp_X(:,3).^2.*temp_X(:,2);
temp_X(:,17)=temp_X(:,1).^0.5;
temp_X(:,18)=temp_X(:,2).^0.5;
temp_X(:,19)=temp_X(:,3).^0.5;

temp_Y(:,1)=temp_PROXY.Qg_MMSCF;

% Calculating inverse matrix to estimate proxy parameters
temp_PARAM=inv(temp_X'*temp_X)*temp_X'*temp_Y;

PROXY_PARAM{i}=temp_PARAM;

% Checking R2 and RMSE

WELLS_R2_Ori(i)=1-(sum((temp_Y-(temp_X*PROXY_PARAM{i})).^2))...
    /(sum((temp_Y-(mean(temp_Y))).^2));
WELLS_RMSE_Ori(i)=sqrt(mean((temp_Y-(temp_X*PROXY_PARAM{i})).^2));

% Clear calculated temporal vectors
clear temp_*

end
% Concatenating
WELLS_NAMES=strcat(WELLS_NAMES, {', '},WELLS_GOR);
% Errors
eje_x = 1:33;

figure
plot(eje_x,WELLS_R2,'*r',eje_x,WELLS_R2_Ori,'.b')
title('R^2')
set(gca,'XTick',[1:33]); 
set(gca,'XTickLabel',WELLS_NAMES);
grid on
legend('Pseudopressure','Original')

figure
plot(eje_x,WELLS_RMSE,'*r',eje_x,WELLS_RMSE_Ori,'.b')
title('RMSE')
set(gca,'XTick',[1:33]); 
set(gca,'XTickLabel',WELLS_NAMES);
grid on
legend('Pseudopressure','Original')

% REFERENCES
% https://nataliaacevedo.com/rmse%E2%80%8A-%E2%80%8A-error-cuadratico-medio/
% https://www.datatechnotes.com/2019/02/regression-model-accuracy-mae-mse-rmse.html
