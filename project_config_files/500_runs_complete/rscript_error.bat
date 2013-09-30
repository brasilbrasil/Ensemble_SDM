D:
cd "D:\Dropbox\code\Ensemble_SDM\project_config_files\500_runs_complete"
::cd to directory where r code is located
::have the number of instances equal the number of CPU cores available
::add plenty of time between initiations to allow for instances to detect total number of running instances

::needs time to copy files in first instance
START "" "C:\Program Files\R\R-3.0.1\bin\x64\Rscript.exe" CMD BATCH --no-save --no-restore --verbose test.r > outputFile.txt 2>&1 /b
exit

START "" "C:\Program Files\R\R-3.0.1\bin\x64\Rscript.exe" --no-save --no-restore --verbose test.r > outputFile.txt 2>&1 /b
RScript --no-save --no-restore --verbose myRfile.R > outputFile.Rout 2>&1

START "" "C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH --vanilla "test.r"

"C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH --vanilla --slave "test.r"

START "" "C:\Program Files\R\R-3.0.1\bin\x64\Rscript.exe" CMD BATCH --vanilla --slave --no-restore --verbose test.r > outputFile.txt 2>&1 /b
"C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH test.r --vanilla --slave --no-restore --verbose test.r > outputFile.txt 2>&1 /b

works!! but no progress
"C:\Program Files\R\R-3.0.1\bin\x64\Rscript.exe" test.r --vanilla --slave --no-restore --verbose test.r > outputFile.txt 2>&1 /b

works!!
"C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH test.r output.txt --vanilla --slave --no-restore --verbose > outputFile.txt 2>&1 /b

Works and cleaner!!
"C:\Program Files\R\R-3.0.1\bin\x64\R.exe" CMD BATCH "test.r" output.txt
