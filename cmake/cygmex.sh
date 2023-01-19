#!/bin/bash
# This is a script of tricks to compile matlab interface with mex and Visual Studio when executed from cygwin

mexargs=""

for arg in $*; do 
	if [[ $arg == *"/"* ]]; then
		if [[ $arg == "-I"* ]]; then
			conv_arg="-I`cygpath -m ${arg:2}`"
		elif [[ $arg == "-L"* ]]; then
			conv_arg="-L`cygpath -m ${arg:2}`"
		else
			conv_arg=`cygpath -m $arg`
		fi
	else
		conv_arg=$arg
	fi	
	if [[ $conv_arg == *".cpp.o" ]]; then
		conv_arg=`echo ${conv_arg%.cpp.o}.obj | sed 's/__\/pfs//g'`				
	fi
	mexargs="$mexargs $conv_arg"
done

#echo "Running command: @MATLAB_MEX_EXE@ $mexargs"

"@MATLAB_MEX_EXE@" $mexargs

 