. "./common.ps1"

try
{
    tar -czf windows.tar.gz -C $INSTALL_DIR .
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

