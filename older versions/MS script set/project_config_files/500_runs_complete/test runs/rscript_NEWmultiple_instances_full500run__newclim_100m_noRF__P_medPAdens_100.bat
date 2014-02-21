D:
cd "D:\Dropbox\code\Ensemble_SDM\project_config_files\500_runs_complete"
::cd to directory where r code is located
::have the number of instances equal the number of CPU cores available
::add plenty of time between initiations to allow for instances to detect total number of running instances

::needs time to copy files in first instance
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_medPAdens_100.r" PA_output_instancePA1.txt /b
timeout 30
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_medPAdens_100.r" PA_output_instancePA2.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_medPAdens_100.r" PA_output_instancePA3.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_medPAdens_100.r" PA_output_instancePA4.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_medPAdens_100.r" PA_output_instancePA5.txt /b

exit
