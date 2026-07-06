@echo off
REM Push this repo to https://github.com/xulinpan/scabm  (double-click or run in cmd/PowerShell)
cd /d "%~dp0"
git init
git config user.name "Xulin Pan"
git config user.email "xulinpanias@gmail.com"
git add -A
git commit -m "SC-ABM submission: manuscript, code, data, NEON application"
git branch -M main
git remote remove origin 2>nul
git remote add origin https://github.com/xulinpan/scabm.git
git push -u origin main
echo.
echo ===== Done. If prompted, sign in to GitHub (use a Personal Access Token as the password). =====
pause
