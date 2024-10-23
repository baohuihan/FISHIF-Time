function poolsize=PoolCreate(NWorkers)
poolobj = gcp('nocreate');% If no pool,  create new one.
if isempty(poolobj)
    poolsize = NWorkers;
    CoreNum=NWorkers; %
    parpool(CoreNum);
else
    poolsize = poolobj.NumWorkers;
    disp('Already initialized'); 
end


