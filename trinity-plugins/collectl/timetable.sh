FILES=`ls *.sum | sort --numeric-sort`

# parse interval of collect
INTERVAL=`cat global.time | grep interval | cut -d ' ' -f 2`
INTLEN=${#INTERVAL}
if [ $INTLEN -gt 1 ]; then
    INTERVAL=`echo ${INTERVAL} | cut -d ':' -f 2`
fi

# start and stop times of chrysalis
FILE=`ls *.Chrysalis.sum`
start_d=`head -n 1 $FILE | cut -d ' ' -f 1-2`
chrysalis_start=`date -d "$start_d" +"%s"`
end_d=`tail -n 1 $FILE | cut -d ' ' -f 1-2`
chrysalis_end=`date -d "$end_d" +"%s"`

chrysalis_counter=0

echo "application runtime[s]" >runtime.csv
for FILE in $FILES; do
if [ -f $FILE ] ; then
  lines=`cat $FILE | wc -l`
  runtime=$(( $lines * $INTERVAL ))
  app_name=`echo $FILE | cut -d '.' -f 3`
# START only for chrysalis
  start_d=`head -n 1 $FILE | cut -d ' ' -f 1-2`
  start_s=`date -d "$start_d" +"%s"`
  end_d=`tail -n 1 $FILE | cut -d ' ' -f 1-2`
  end_s=`date -d "$end_d" +"%s"`
  # if tool inside of chrysalis - add to counter 
  if [ $start_s -ge $chrysalis_start -a  $end_s -le $chrysalis_end ]; then
    if [ "$app_name" != "Chrysalis" ]; then
      echo "$app_name inside Chrysalis"
      chrysalis_counter=$(($chrysalis_counter + $runtime))
    fi
  fi
  # count only exclusive time spent in chrysalis
  if [ "$app_name" = "Chrysalis" ]; then
    runtime=$(($runtime - $chrysalis_counter))
  fi
# END only for chrysalis
  echo "$app_name $runtime" >>runtime.csv
fi
done

