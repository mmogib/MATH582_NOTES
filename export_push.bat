@echo off
setlocal

rem %~1 = %1 with surrounding quotes removed
set "MESSAGE=%~1"

julia --project=. src\export.jl || goto :fail
git add .
git commit -m "%MESSAGE%" || goto :fail
git push || goto :fail
goto :eof

:fail
echo Failed with error %errorlevel%.
exit /b %errorlevel%
