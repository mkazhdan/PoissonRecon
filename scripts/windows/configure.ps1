. "./common.ps1"

try
{
    powershell -nologo -command "cmake -B $BUILD_DIR -S $PROJECT_ROOT $CMAKE_FLAGS"
}
catch [System.SystemException]
{
    Write-Output "Error occured:"
    Write-Output $_
}
finally
{
    cd $SCRIPTS_DIR
}

