function jobLoc = launch_parpool(nWorkers, A, a)
%% launch parpool

jobLoc = sprintf('~/scratch/matlabJobs/%d_%d', A, a)


if exist(jobLoc, 'dir')
    rmdir(jobLoc, 's');
end
mkdir(jobLoc);


myCluster = parcluster;
myCluster.JobStorageLocation = jobLoc;
myCluster.NumThreads = 2*nWorkers;

parpool(myCluster, nWorkers)