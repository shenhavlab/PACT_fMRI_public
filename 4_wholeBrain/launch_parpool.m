function launch_parpool(nWorkers)
%% launch parpool

myCluster = parcluster;
myCluster.JobStorageLocation = '~/scratch/matlabJobs';
myCluster.NumThreads = 2*nWorkers;

parpool(myCluster, nWorkers)