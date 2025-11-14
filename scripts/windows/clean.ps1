. "./common.ps1"

try
{
    rm -r $BUILD_DIR
    rm -r $INSTALL_DIR
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

