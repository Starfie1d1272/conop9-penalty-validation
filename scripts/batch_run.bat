@echo off
chcp 65001 >nul

rem ---- repo-relative paths（从脚本位置推导仓库根）----
set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "REPO_DIR=%%~fI"
set "CONOP_DIR=%REPO_DIR%\CONOP-run"
set "RESULTS_DIR=%REPO_DIR%\results"
set REPEAT=3

rem ---- external executable policy ----
rem 仓库不提供 CONOP9 executable。默认按文件名提示，用户可先 set CONOP_EXE=本地路径
if not defined CONOP_EXE set "CONOP_EXE=%CONOP_DIR%\CONOP64ver8p621.exe"
if not exist "%CONOP_EXE%" (
    echo CONOP9 executable is not distributed with this repository.
    echo Set CONOP_EXE to the path of your local CONOP9 executable, e.g.:
    echo   set CONOP_EXE=C:\path\to\CONOP64ver8p621.exe
    echo 仓库内未找到 executable；请先设置 CONOP_EXE 再运行本脚本。
    exit /b 1
)

cd /d "%CONOP_DIR%"

echo ============================================
echo   CONOP 批量参数扫描脚本（每组 %REPEAT% 次）
echo ============================================
echo   executable: %CONOP_EXE%
echo.

call :run_one baseline      0.980  250  600
call :run_one ratio_099    0.990  250  600
call :run_one ratio_095    0.950  250  600
call :run_one temp_500     0.980  500  600
call :run_one temp_100     0.980  100  600
call :run_one steps_1200   0.980  250 1200
call :run_one steps_0300   0.980  250  300

echo.
echo ===== 全部完成！结果在 %RESULTS_DIR%\ =====
echo.
echo git add -A
echo git commit -m "批量扫描结果（每组重复3次）"
echo git push
echo.
pause
exit /b

:run_one
set TAG=%1
set NEW_RATIO=%2
set NEW_TEMP=%3
set NEW_STEPS=%4

for /l %%i in (1,1,%REPEAT%) do (
    echo ----------------------------------------
    echo [%TAG%] 第 %%i 次: RATIO=%NEW_RATIO%  TEMP=%NEW_TEMP%  STEPS=%NEW_STEPS%
    echo.

    powershell -Command "(Get-Content conop9.cfg) -replace 'RATIO=[0-9.]+', 'RATIO=%NEW_RATIO%' | Set-Content conop9.cfg" 2>nul
    powershell -Command "(Get-Content conop9.cfg) -replace 'STARTEMP=[0-9.]+', 'STARTEMP=%NEW_TEMP%' | Set-Content conop9.cfg" 2>nul
    powershell -Command "(Get-Content conop9.cfg) -replace 'STEPS=[0-9]+', 'STEPS=%NEW_STEPS%' | Set-Content conop9.cfg" 2>nul

    echo 配置已更新，请用 %CONOP_EXE% 指向的本地 executable 运行
    echo 跑完后关掉窗口，按任意键继续...
    echo.
    pause

    :: 结果统一存到仓库相对 results 目录
    mkdir "%RESULTS_DIR%\%TAG%\run_%%i" 2>nul

    copy /Y trajectory.txt  "%RESULTS_DIR%\%TAG%\run_%%i\trajectory.txt"  >nul
    copy /Y bestsoln.dat    "%RESULTS_DIR%\%TAG%\run_%%i\bestsoln.dat"    >nul
    copy /Y outmain.txt     "%RESULTS_DIR%\%TAG%\run_%%i\outmain.txt"     >nul
    copy /Y runlog.txt      "%RESULTS_DIR%\%TAG%\run_%%i\runlog.txt"      >nul
    copy /Y conop9.cfg      "%RESULTS_DIR%\%TAG%\run_%%i\conop9.cfg"      >nul
    copy /Y soln.dat        "%RESULTS_DIR%\%TAG%\run_%%i\soln.dat"        >nul
    copy /Y ordr.dat        "%RESULTS_DIR%\%TAG%\run_%%i\ordr.dat"        >nul

    if exist "%RESULTS_DIR%\%TAG%\run_%%i\trajectory.txt" (
        echo [%TAG%] 第 %%i 次 - 保存成功
    ) else (
        echo [%TAG%] 第 %%i 次 - 保存失败
        echo 当前目录：
        dir
    )
)

echo ----------------------------------------
exit /b
