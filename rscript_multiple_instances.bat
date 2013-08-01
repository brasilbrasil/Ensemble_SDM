cd "Y:\PICCC_analysis\code\Ensemble_SDM"
::cd to directory where r code is located
::have the number of instances equal the number of CPU cores available
::add plenty of time between initiations to allow for instances to detect total number of running instances
START "" "C:\Program Files\R\R-3.0.0\bin\x64\Rscript.exe" "1_BM2_FB_SDM_fitting_w_cropping4.r" /b
timeout 60
START "" "C:\Program Files\R\R-3.0.0\bin\x64\Rscript.exe" "1_BM2_FB_SDM_fitting_w_cropping4.r" /b
timeout 10
START "" "C:\Program Files\R\R-3.0.0\bin\x64\Rscript.exe" "1_BM2_FB_SDM_fitting_w_cropping4.r" /b
timeout 10
START "" "C:\Program Files\R\R-3.0.0\bin\x64\Rscript.exe" "1_BM2_FB_SDM_fitting_w_cropping4.r" /b
exit
