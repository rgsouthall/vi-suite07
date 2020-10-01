set NGSCXX_DIR=%~dp0
 call "C:/Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/VC/Auxiliary/Build/vcvarsall.bat" amd64
link /DLL %*  /LIBPATH:"%NGSCXX_DIR%/../lib" nglib.lib ngcore.lib libngsolve.lib
