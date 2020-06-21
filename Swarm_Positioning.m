
%%
clear all;
va=0.5;

for va=0.01:0.01:1
    
%Specifying Parameters
format long
nx=200;                %Number of steps in space(x)
nt=60000;               %Number of time steps 
dt=0.001; 
%Width of each time step
dx=1/(nx-1);          %Width of space step
x=0:dx:1;             %Range of x (0,2) and specifying the grid points
CB=zeros(nx,1);       %Preallocating u
CR=zeros(nx,1);
CA=zeros(nx,1); 
CT=zeros(nx,1);  
FA=zeros(nx,1);
FR=zeros(nx,1);
CM=zeros(nx,1);
CP=zeros(nx,1);
SecretRate=1;
JSD_timeline=[];
UnitedV1Rate=1;
UnitedV2Rate=1;

DiffA=0.3;            %Diffusion coefficient/viscosity
DiffR=0.05;
DiffB=1;

KA=0.1;               %Mortality
KR=0.1;               %Mortality

r=dt/(dx*dx);

%Set='Uniform';
Set='Left';
vr=0;%0.0001;             %fa

ka=1;                 %fa
kr=1;                 %fr


UL=0;                 %Left Dirichlet B.C
UR=0;                 %Right Dirichlet B.C
UnL=0;                %Left Neumann B.C (du/dn=UnL) 
UnR=0;                %Right Neumann B.C (du/dn=UnR) 

%%%%%%%%%%%%%%%%%%%%%%%%%


CT=(1/2)*normpdf(x, 0.8, 0.05)'+(1/2)*normpdf(x, 0.2, 0.05)';   
% 

% CB=normpdf(x, 0.1, 0.05)';   

for i=1:nx+1
     if (i==15)
        CB(i)=nx;
    end
end

% for i=1:nx
%         CB(i)=1;   
% end


 %CB=normpdf(x, 0.1, 10)';
%CT=(1/2)*(1/(sqrt(2*pi)*0.01))*gaussmf(x, [0.01 0.8])'+(1/2)*(1/(sqrt(2*pi)*0.01))*gaussmf(x, [0.01 0.3])';        
%CB=(1/(sqrt(2*pi)*0.1))*gaussmf(x, [0.1 0.1])';

% for i=1:nx
%     CR(i)=0;
%     CA(i)=0;
%     CB(i)=2*exp(-200*(x(i)-0.5).^2);
%  %   CB(i)=1*exp(-50*(x(i)-1).^2);
% 
%  %%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$               
% 
% %     if ((0.300<=x(i))&&(x(i)<=0.500))
% %         CT(i)=2;
% %       else if ((0.600<=x(i))&&(x(i)<=0.800))
% %        CT(i)=2;
% %           else
% %             if ((1.0<=x(i))&&(x(i)<=1.100))
% %              CT(i)=2;
% %             else
% %                 CT(i)=0;
% %             end
% %           end
% %     end
% % 
% %     if ((0.200<=x(i))&&(x(i)<=0.400))
% %        CT(i)=2;
% %       else 
% %        CT(i)=0;
% %     end
% %     
%     %if ((0.300<=x(i))&&(x(i)<=0.500))CT(i)=(1/sqrt(2*pi))*exp(-((x(i)-0.4).^2)/(2));%(1/sqrt(2*pi))*exp(-((x(i)-0.7).^2)/(2*0.04));
%     
% 
%      % else if ((1.7<=x(i))&&(x(i)<=1.900))
%       %    else
%         %       CT(i)=0;
%          % end
%  end    
%  %%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$               

 

%% fixed part of coefficient matrix of the first equation     No mistake
L=-DiffB*r;  
M=2*DiffB*r;
N=-DiffB*r;
E=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
E=L*E;
Q=sparse(1:nx-3,2:nx-2,1,nx-2,nx-2);
Q=N*Q;
EQ=E+Q;

Y=speye(nx-2);
Y=Y+M*Y;
EYQ=EQ+Y;

EYQ(1,1)=EYQ(1,1)+L;
EYQ(nx-2,nx-2)=EYQ(nx-2,nx-2)+L;

%fixed part for the first equation

%% fixed part of coefficient matrix of the second equation   No mistake
L1=-DiffR*r;
M1=2*DiffR*r+KR*dt;
N1=L1;
E1=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
E1=L1*E1;
Q1=sparse(1:nx-3,2:nx-2,1,nx-2,nx-2);
Q1=N1*Q1;
EQ1=E1+Q1;
Y1=speye(nx-2);
Y1=Y1+M1*Y1;
EYQ1=EQ1+Y1;
EYQ1(1,1)=EYQ1(1,1)+L1;
EYQ1(nx-2,nx-2)=EYQ1(nx-2,nx-2)+L1;

% for the second equation, there is no varied part of its coefficient
% matrix

%% fixed part of coefficient matrix of the third equation
L2=-DiffA*r; 
M2=2*DiffA*r+KA*dt;
N2=L2;
E2=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
E2=L2*E2;
Q2=sparse(1:nx-3,2:nx-2,1,nx-2,nx-2);
Q2=N2*Q2;
EQ2=E2+Q2;
Y2=speye(nx-2);
Y2=Y2+M2*Y2;
EYQ2=EQ2+Y2;
EYQ2(1,1)=EYQ2(1,1)+L2;
EYQ2(nx-2,nx-2)=EYQ2(nx-2,nx-2)+L2;

% for the second equation, there is no varied part of its coefficient
% matrix

itn=zeros();

     JSD_timeline=[];
    %% interior iterator
     for it=1:nt

    %      if it==60
    %          UnitedV1Rate=1;
    %      end

        %Performance Record
        CAtemp=CA;
        CBtemp=CB;
        CRtemp=CR;

    %        M=(CB+CT)/2;
    %        JSD=real((sum(CB.*log2(CB./M)/(nx+1))+sum(CT.*log2(CT./M)/nx))/2);
    %        JSD_timeline=[JSD_timeline; JSD];




           M=(CB+CT)/2;
           JSD=real((sum(CB.*log2(CB./M)/(nx+1))+sum(CT.*log2(CT./M)/nx))/2);
           JSD_timeline=[JSD_timeline; JSD];





    %        if it==5001
    %            pause;
    %         filename0=['C:\Users\jian\Documents\MATLAB\','Set=', Set, 'JSD, Va=',num2str(va),',Vr=', num2str(vr),'.csv']; 
    %         xlswrite(filename0, JSD_timeline);
    %         load train;
    %         sound(y,Fs);
    %         return;
    %        end


    % 
    %         filename2=['C:\Users\jian\Documents\MATLAB\TempTimeSeries\', 'Set=', Set, 'Va=',num2str(va),',Vr=',num2str(vr),',CA_t=',num2str(it*dt*10),'.csv'];
    %         csvwrite(filename2, CA);
    % %         filename3=['C:\Users\jian\Documents\MATLAB\TempTimeSeries\', 'Set=', Set, 'Va=',num2str(va),',Vr=',num2str(vr),',CR_t=',num2str(it*dt*10),'.csv'];
    % %         csvwrite(filename3, CR);
    %         filename4=['C:\Users\jian\Documents\MATLAB\TempTimeSeries\', 'Set=', Set, 'Va=',num2str(va),',Vr=',num2str(vr),',CB_t=',num2str(it*dt*10),'.csv'];
    %         csvwrite(filename4, CB);


      %refreshdata(h)
    % refreshdata(j)


        %Uncomment as necessary
        %-------------------

            CAt=CAtemp; CBt=CBtemp; CRt=CRtemp;  CTt=CT;
            CAt(1)=[];CAt(end)=[];  CBt(1)=[];CBt(end)=[];  CRt(1)=[];CRt(end)=[];   CTt(1)=[]; CTt(end)=[];
         %Implicit solution for the second equation      No mistake
    %         RRightPart=SecretRate*CBt*dt+CRt;
    %         CRt=EYQ1\RRightPart;
    %         CRtemp=[CRt(1);CRt;CRt(end)]; %Neumann
    %         CR=CRtemp;
    %         
            Temp=(dt*UnitedV2Rate*(CBt.*(1./(CTt+kr))));    %./((UnitedK1Rate+CBt).*(UnitedK2Rate+CTt));
            RRightPart=Temp+CRt;
            CRt=EYQ1\RRightPart;
            CRtemp=[CRt(2);CRt;CRt(end-1)]; %Neuman
            CR=CRtemp;
    % 

         %Implicit solution for the third equation      No mistake
            Temp=(dt*UnitedV1Rate*(CBt.*(CTt./(CTt+ka))));         %./((UnitedK1Rate+CBt).*(UnitedK2Rate+CTt));
            ARightPart=Temp+CAt;
            CAt=EYQ2\ARightPart;
            CAtemp=[CAt(2);CAt;CAt(end-1)]; %Neuman
            CA=CAtemp;
        %Implicit solution for the second and the third equation finished


        %Implicit solution for the first equation 

            %varied part of the first equation start   No mistake
            VYQ=zeros(nx-2,nx-2);% declaration

            for i=2:nx-1
                FA(i)=(CA(i+1)-CA(i-1))/(2*dx);
                %FA(i)=(FA(i)+0.01)*(va*CA(i)/(ka+CA(i)));
                FR(i)=(CR(i+1)-CR(i-1))/(2*dx);
                %FR(i)=(FR(i)+0.01)*(vr*CR(i)/(kr+CR(i)));
            end

        F=(va*FA-vr*FR);
        %    for i=2:nx-1
          %      if CA(i)>CR(i)

     %%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$               
                 %    F(i)=0.06*(FA(i)-FR(i));
                 %        F(i)=0.008*(FA(i)-FR(i));

               % else 
                  %     F(i)=1.0*(FA(i)-FR(i));
     %%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$               
              %  end
           % end

    %         for i=2:nx-1
    %             if CR(i)>CA(i)
    %                 F(i)=(FR(i)-FA(i));
    %             else if CR(i)<=CA(i)
    %                 
    %                 F(i)=-(FR(i)-FA(i));
    %                 end
    %             end
    %         end


            for i=1:nx-2
                for j=1:nx-2
                    if j==i+1
                        VYQ(i,j)=VYQ(i,j)+dx*F(i+2)*r/2;
    %                      VYQ(i,j)=VYQ(i,j)+dx*0.3*r/2;
                    end
                    if j==i-1
                        VYQ(i,j)=VYQ(i,j)-dx*F(i)*r/2;
    %                     VYQ(i,j)=VYQ(i,j)-dx*3*r/2;
                    end
                end
            end



            %varied part of the first equation finished   No mistake

            %combine to form the coeffcient matrix  No mistake
            EYQ=EYQ+VYQ;
            %EYQ(1,1)=EYQ(1,1)-dx*F(2)*r/2;
            %EYQ(nx-2,nx-2)=EYQ(nx-2,nx-2)+dx*F(nx-1)*r/2;
            %EYQ(1,2)=EYQ(1,2)-dx*F(2)*r/2;
            %EYQ(nx-2,nx-3)=EYQ(nx-2,nx-3)+dx*F(nx-1)*r/2;

            %This is the coefficient matrix for the first partial differential equation
            CBt=EYQ\CBt;
            CBtemp=[CBt(2);CBt;CBt(end-1)]; %Neumann
            CB=CBtemp;


     end
     
     
      [pks,locs] = findpeaks(JSD_timeline(40001:60000));
      % meanCycleA = mean(diff(locs));
      [pks1,locs1] = findpeaks(-JSD_timeline(40001:60000));
       pks1=-pks1; locs1=-locs1;     
      % meanCycleB = mean(diff(locs1));
      
      if isempty(pks)
          pks=JSD_timeline(50000);
          pks1=JSD_timeline(50000);
          
      end
%       
%       
%      
%       
%        
%        filename=['D:\Simu_Yao\', 'pks____Set=', Set, ',Va= ',num2str(va),',Vr= ', num2str(vr),'.txt']; 
%        fileID = fopen(filename,'w');
%        fprintf(fileID,'%f',pks);
%        fclose(fileID);
%        
%        filename=['D:\Simu_Yao\', 'pks1____Set=', Set, ',Va= ',num2str(va),',Vr= ', num2str(vr),'.txt']; 
%        fileID = fopen(filename,'w');
%        fprintf(fileID,'%f',pks1);
%        fclose(fileID);
       
%        filename=['D:\Simu_Yao\', 'locs____Set=', Set, ',Va= ',num2str(va),',Vr= ', num2str(vr),'.txt']; 
%        fileID = fopen(filename,'w');
%        fprintf(fileID,'%f',locs);
%        fclose(fileID);
% 
%        filename=['D:\Simu_Yao\', 'locs1____Set=', Set, ',Va= ',num2str(va),',Vr= ', num2str(vr),'.txt']; 
%        fileID = fopen(filename,'w');
%        fprintf(fileID,'%f',locs1);
%        fclose(fileID);
       
%        filename=['D:\Simu_Yao\', 'JSD_timeline____Set=', Set, ',Va= ',num2str(va),',Vr= ', num2str(vr),'.txt']; 
%        fileID = fopen(filename,'w');
%        fprintf(fileID,'%f',JSD_timeline);
%        fclose(fileID);
%      
       sizepks=length(pks);
       sizepks1=length(pks1);
       vrarray=va*ones(sizepks,1);
       vrarray1=va*ones(sizepks1,1);
       y=plot(vrarray,pks,'ob');
       hold on;
       y=plot(vrarray1,pks1,'ob');
       hold on;
%        
      
end
                      


axis([0 1 0 2])
drawnow;



        figure(2);
    %     subplot(121)
        subplot(121)
       % set (gca,'position',[0.1,0.1,0.9,0.9] );
        h=plot(x,CR,'g','LineWidth',2);       %plotting the concentration profile
        hold on;
        h=plot(x,CA,'r','LineWidth',2);
        hold on;
        h=plot(x,CB,'k','LineWidth',2);
        hold on;
        h=plot(x,CT,'b','LineWidth',2);
        hold on;
        %h=imresize(h,[256,256]);
        axis square
        axis([0 1 0 15])

           title({['1-D simulation'];['Time(\itt) = ',num2str(dt*it*10)];  ['Initial condition = ', Set]; ['V a = ',num2str(va),',  V r = ',num2str(vr)];['B=',num2str(sum(CB/nx))]})   %;['Norm = ',num2str(L2norm(it))] 

            % ['B=',num2str(sum(CB/nx))]; 
        xlabel('Spatial coordinate (x) \rightarrow')
        ylabel('Concentration profile \rightarrow')
     %   drawnow; 

    %     filename1=['C:\Users\jian\Documents\MATLAB\', 'Set=', Set, 'Va= ',num2str(va),',Vr= ', num2str(vr),'.gif']; 
    %     frame = getframe(1);
    %       im = frame2im(frame);
    %       [imind,cm] = rgb2ind(im,256);
    %       if it == 1
    %          imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
    %       else
    %          imwrite(imind,cm,filename1,'gif','WriteMode','append');
    %       end
    % 



    %      hold off;
    %     drawnow;
       %figure(2);
        dp=0.01;
        p=0:dp:(it-1)/100;
        subplot(122)
        j=plot(p, JSD_timeline,'LineWidth',1.5);
        axis square
            axis([400 600 0 1])
           title(['Real-time JSD']);
            %title(['\itJSD= ',num2str(JSD)]);
        xlabel('Time (t) \rightarrow')
        ylabel('JSD value')
        drawnow;


         refreshdata(h);
         refreshdata(j);
% 



