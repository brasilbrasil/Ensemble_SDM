D:
cd "D:\Dropbox\code\Ensemble_SDM\project_config_files\500_runs_complete"
::cd to directory where r code is located
::have the number of instances equal the number of CPU cores available
::add plenty of time between initiations to allow for instances to detect total number of running instances

::needs time to copy files in first instance
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance1laptop.txt /b
timeout 30
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance2laptop.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance3laptop.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance4laptop.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance5laptop.txt /b
timeout 10
START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "model_runs_with_old_code_noAbs.r" final_output_instance6laptop.txt /b

exit
