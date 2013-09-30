D:
cd "D:\Dropbox\code\Ensemble_SDM\project_config_files\500_runs_complete"
::cd to directory where r code is located
::have the number of instances equal the number of CPU cores available
::add plenty of time between initiations to allow for instances to detect total number of running instances

::needs time to copy files in first instance
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance1.txt /b
timeout 30
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance2.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance3.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance4.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance5.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance6.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance7.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance8.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance9.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance10.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance11.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "new_code_full500run__newclim_100m_noRF__P_A_medPAdens_1.r" output_instance12.txt /b

exit
