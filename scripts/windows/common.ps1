# Exit on error:
$ErrorActionPreference = "Stop" 

Set-Variable -Name PLATFORM -Value "windows"
Set-Variable -Name SCRIPTS_DIR -Value (Get-Location)
Set-Variable -Name PROJECT_ROOT -Value $SCRIPTS_DIR\..\..\
Set-Variable -Name INSTALL_DIR -Value $PROJECT_ROOT\install\$PLATFORM
Set-Variable -Name BUILD_DIR -Value $PROJECT_ROOT\build\$PLATFORM
# FIXME
Set-Variable -Name CPU_COUNT -Value 8
Set-Variable -Name BUILD_TYPE -Value Release

# Settle on specific toolchain version:
$CMAKE_FLAGS = " -G 'Visual Studio 17 2022'"
$CMAKE_FLAGS = " -T 'v143'"
$CMAKE_FLAGS += " -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR"
$CMAKE_FLAGS += " -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL"

New-Item $INSTALL_DIR -ItemType Directory -Force
New-Item $BUILD_DIR -ItemType Directory -Force

