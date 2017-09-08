if ( ! $?LD_LIBRARY_PATH ) then
   echo "Setting LD_LIBRARY_PATH"
   setenv LD_LIBRARY_PATH $PWD/source/lib
else
   echo "Updating LD_LIBRARY_PATH"
   setenv LD_LIBRARY_PATH $PWD/source/lib:$LD_LIBRARY_PATH    
endif
echo $LD_LIBRARY_PATH
