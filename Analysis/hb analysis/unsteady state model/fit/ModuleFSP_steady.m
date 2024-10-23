function [lld,varargout] = ModuleFSP_steady(Model,Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin)
%%% A program to save and show the fitting parameter
switch Model
    case 'TwoState'
        [lld,varargout]=FSP_FreeSet_Steady(Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin);
    case 'Ton'
        [lld,varargout]=FSP_FreeSetTon_Steady(Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin);
end