function [lld,varargout] = ModuleFSP(Model,Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin)
%%% A program to save and show the fitting parameter
switch Model
    case 'TwoState'
        [lld,varargout]=FSP_FreeSet(Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin);
    case 'Ton'
        [lld,varargout]=FSP_FreeSetTon(Kanneal,ndata,p_ini,dna_label,hold_time,Kall,AnnealLogic,Tstay,DtdInherit,varargin);
end