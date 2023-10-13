set i=13
set f=298

while ( $i < $f )
 
  echo "Now Working on: " $i "   File"
 
  sh run_MDC2reco.sh FileList/${i}.list
  
  @ i = $i + 1 

end
