. "./common.ps1"

try
{
    powershell -nologo -command "cmake --build $BUILD_DIR --config $BUILD_TYPE --target install -j $CPU_COUNT"
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

