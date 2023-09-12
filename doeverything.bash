source ../trapenvvar.bash
cp ../../deduplicate.py . 
python deduplicate.py --dataset 4 --pullsql --getoutput
cp ../../slidingetawindow.py . 
python slidingetawindow.py "${TKP_DBNAME}deduped.csv"
cp ../../commensal2trig*csv .
cp ../../makelcfromruncat.py . 
python makelcfromruncat.py *sorted*csv
