function [lld,varargout] = FSP_FreeSetPon2(Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,Change,DtdInherit,varargin)
%% A program to calculate the likelihood of 3 state model
% All mrna number can be OFF
%% Parameter setting and initialization
Kall(AnnealLogic)=Kanneal;
k_all=Kall;
T = dna_label/sum(dna_label);
if isempty(varargin) || isempty(varargin{1})
    Nmax = 100;
    TXmax = 100;
else
    [Nmax,Nmax2] = varargin{1};
end
Nstate=3;
Nbin=1; %1 mRNA1
TimeBin=10;
% Nbin2=1;%1 mRNA2
Nbint=Tstay*TimeBin; %total time 
Nbint_pre1 = floor(Nbint*T(1)); % Number of time bin 
Nbint_dur1 = floor(Nbint*T(2));    
Nbint_post1 = floor(Nbint*T(3));
Nbint_post1 = Nbint_post1 + floor(Nbint*hold_time);
%dt2 = t2-Nbint2/Nbint;
%dt3 = t3-Nbint3/Nbint;
Nr = Nbint_dur1/Nbin;
dt = 1/TimeBin;                %%% time step 0.1
dRNA = 1/Nbin;               %%% length step
% dRNA2 = 1/Nbin2;
nRNA = 0:dRNA:Nmax;          %%% number of complete RNA molecules
% nRNA2 = 0:dRNA2:Nmax2;
Nx = length(nRNA);           %%% max mRNA step number 
% Nx2 = length(nRNA2);
TimeEvo=0;
%% k_tranition
k_intr=k_all(end-2)*dt;% useless now
k=k_all(1:9);
inten_intr=k_all(end-1);
GeneActiveR=k_all(end-2);
ndata=ndata/inten_intr;
k_tra = [-(k(1)+k(2)),k(3),k(5);k(1),-(k(3)+k(4)),k(6);k(2),k(4),-(k(5)+k(6))];
k_tra_p=k_tra*dt+diag(ones(1,3));
km1 = diag(k(7:9));
%% k_matirx
k_self1 = k_tra*dt+(diag(ones(1,Nstate))-km1*dt);
k_smat1 = diag_matirx_1(k_self1,Nx,Nx);
k_smat1_use = k_smat1;
%% initialization
% p_ini = [sum(k_tra(1,:))-k_tra(1,1),sum(k_tra(2,:))-k_tra(2,2),sum(k_tra(3,:))-k_tra(3,3),sum(k_tra(4,:))-k_tra(4,4)]'/sum(k(1:12));
if  isempty(Change)
    Vary=GeneActiveR-sum(p_ini);
    if Vary>=0
        p_ini(1)=p_ini(1)+Vary;
   else
        p_ini=p_ini+p_ini*Vary/sum(p_ini);
   end
end

Change=[Change,[inf;GeneActiveR]];
Change(1,:)=Change(1,:)*TimeBin;
% if  ~isempty(Change)
%     p_ini=p_ini*Change(2,1);
% else
%     p_ini=p_ini*GeneActiveR;
% end
ChangeVary=Change;
if size(ChangeVary,2)>1
    ChangeVary(2,2:end)=ChangeVary(2,2:end)-ChangeVary(2,1:end-1);
end

Error=0;
mRNAadd=zeros(1,Nx);
if isempty(DtdInherit)
    P_dtb_1 = zeros(Nstate*Nx,1);
    P_dtb_1(1:Nstate,1)=p_ini;
else
    P_dtb_1=DtdInherit;
    Vary=GeneActiveR-sum(sum(P_dtb_1));
   if Vary>=0
        P_dtb_1(1)=P_dtb_1(1)+Vary;
   else
       Pdel=P_dtb_1*Vary/sum(P_dtb_1);
       P_dtb_all=reshape(Pdel,[Nstate,Nx]);
        mRNAadd =mRNAadd-sum(P_dtb_all,1);
        P_dtb_1=P_dtb_1+P_dtb_1*Vary/sum(P_dtb_1);
   end
   if sum(P_dtb_1<0)>0
       Error=inf;
   end
end

%% post probe
rna_signal1=Nbin;
k_smat1_use = k_smat1;
[x_other,y_other]=size(k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1))));
k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))=k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))+diag_matirx_1(km1*dt,x_other/Nstate,y_other/Nstate);
for t = Nbint_post1:-1:1
    P_dtb_1=k_smat1_use*P_dtb_1;
    TimeEvo=TimeEvo+1;
    if ~isempty(Change)
        if TimeEvo>Change(1,1)
           Vary=ChangeVary(2,2);
           if Vary>=0
                P_dtb_1(1)=P_dtb_1(1)+ChangeVary(2,2);
           else
               Pdel=P_dtb_1*Vary/sum(P_dtb_1);
               P_dtb_all=reshape(Pdel,[Nstate,Nx]);
                mRNAadd =mRNAadd-sum(P_dtb_all,1);
                P_dtb_1=P_dtb_1+P_dtb_1*Vary/sum(P_dtb_1);
           end
           if sum(P_dtb_1<0)>0
               Error=inf;
           end
           Change(:,1)=[];
           ChangeVary(:,1)=[];
        end
    end
end
rna_signal1_bf=rna_signal1;
%% during probe
k_smat1_use = k_smat1;
[x_other,y_other]=size(k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1))));
k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))=k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))+diag_matirx_1(km1*dt,x_other/Nstate,y_other/Nstate);
for t = Nbint_dur1:-1:1
%     rna_signal1=round((t-1)/Nr);
    rna_signal1=1;
    if rna_signal1~=rna_signal1_bf
        k_smat1_use = k_smat1;
        [x_other,y_other]=size(k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1))));
        k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))=k_smat1_use((Nstate*(rna_signal1)+1):end,1:end-(Nstate*(rna_signal1)))+diag_matirx_1(km1*dt,x_other/Nstate,y_other/Nstate);
        rna_signal1_bf=rna_signal1;
    end
    P_dtb_1=k_smat1_use*P_dtb_1;
    TimeEvo=TimeEvo+1;
    if ~isempty(Change)
        if TimeEvo>Change(1,1)
           Vary=ChangeVary(2,2);
           if Vary>=0
                P_dtb_1(1)=P_dtb_1(1)+ChangeVary(2,2);
           else
                Pdel=-1*P_dtb_1*Vary/sum(P_dtb_1);
               P_dtb_all=reshape(Pdel,[Nstate,Nx]);
                mRNAadd =mRNAadd-sum(P_dtb_all,1);
                P_dtb_1=P_dtb_1+P_dtb_1*Vary/sum(P_dtb_1);
           end
           if sum(P_dtb_1<0)>0
               Error=inf;
           end
           Change(:,1)=[];
           ChangeVary(:,1)=[];
        end
    end
end 

if size(Change,2)>1
     Vary=sum(ChangeVary(2,2:end),2);
       if Vary>=0
            P_dtb_1(1)=P_dtb_1(1)+Vary;
       else
            P_dtb_1(1:3)=P_dtb_1(1:3)+P_dtb_1(1:3)*Vary/sum(P_dtb_1(1:3));
       end
       if sum(P_dtb_1<0)>0
           Error=inf;
       end
       Change(:,1)=[];
       ChangeVary(:,1)=[];
end

if sum(P_dtb_1<0)>0
   Error=inf;
end
P_dtb_all=reshape(P_dtb_1,[Nstate,Nx]);
P_dtb_rna = sum(P_dtb_all);
p_ini=sum(P_dtb_all,2);
P_dtb_nx=reshape(P_dtb_rna(2:end),[Nbin,Nmax]);
P_dtb_rna=[P_dtb_rna(1),sum(P_dtb_nx,1)];
%% pre probe
% no signal
%% double copy
alp=k_all(end);
% P_dtb_rna_single=(1-alp)*[0,P_dtb_rna(2:end)];
% P_dtb_rna_nz=[0,P_dtb_rna(2:end)];P_dtb_rna_nz=P_dtb_rna_nz/sum(P_dtb_rna_nz);
% P_dtb_double=conv(P_dtb_rna_nz,P_dtb_rna_nz);
% P_dtb_rna_double=alp*P_dtb_double;
% P_dtb_rna=[P_dtb_rna(1),zeros(1,length(P_dtb_rna)-1)]+(1-P_dtb_rna(1)).*(P_dtb_rna_single+P_dtb_rna_double(1:length(P_dtb_rna_single)));
P_dtb_rna=P_dtb_rna+mRNAadd;
P_dtb_rna_cov=P_dtb_rna/sum(P_dtb_rna);
p_conv=alp*sum(P_dtb_rna)*conv(P_dtb_rna_cov,P_dtb_rna_cov);
if isnan(p_conv(1))
    P_dtb_rna=(1-alp)*P_dtb_rna;
else  
    P_dtb_rna=(1-alp)*P_dtb_rna+p_conv(1:length(P_dtb_rna));
end
%% cut out the range that experiment is not dependable
P_dtb_rna = P_dtb_rna(1:TXmax+1);
P_dtb_rna(1)=P_dtb_rna(1)+1-sum(P_dtb_rna);
% P_dtb_rna=(GeneActiveR).*P_dtb_rna;P_dtb_rna(1)=P_dtb_rna(1)+1-GeneActiveR;

P_dtb_rna_use=P_dtb_rna(1:71);
P_dtb_rna_use(1)=sum(P_dtb_rna_use(1:4));
P_dtb_rna_use(2:4)=0;
ndata(ndata >= 70.5)=0;
ndata(ndata < 3.5)=0;
%% likelihood
p_tot_bin_log=log(P_dtb_rna_use);
p_tot_bin_log(P_dtb_rna_use<0)=-inf;
lld=-sum(p_tot_bin_log(round(ndata)+1))+Error;%min
% Nbin_data=0:1:100;
% nf = histc(reshape(ndata,1,length(ndata)),Nbin_data);
% nf=nf'/sum(nf);
% lld=-sum(1./abs(P_dtb_rna(1:61)-nf(1:61)'));
varargout={P_dtb_rna,P_dtb_all,p_ini,P_dtb_1};