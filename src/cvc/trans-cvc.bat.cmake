@echo off

set PATH=${BOOST_LOCATION};%PATH%
${BINARIES_LOCATION}\bin\Debug\trans-cvc.exe %*

