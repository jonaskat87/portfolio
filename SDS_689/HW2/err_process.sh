# bash script to process the outputs of 1.py, 2.py, and 3.py, 
# then input into plt_rel_errors.py
# these have been piped into scripts of the forms
# 1_run*.txt, 2_run*.txt, and 3_run*.txt, respectively
cd /c/Users/jonas/Documents/SciML_class/HW2 # change to correct local folder

rm tmp0.txt
rm tmp.txt
# get widths tested (assume same for all cases), save in tmp.txt
grep -i 'L2' 1_run_1.txt | awk '{print $6}' | sed s/://g > tmp0.txt
# load samples for each sample
cp tmp0.txt tmp.txt
# loop through all errors 
for file in 3_run_*.txt # change number depending on the problem (1, 2, or 3)
do 
    # add column with new L2 relative error samples
	paste tmp0.txt <(grep -i 'L2' $file | awk '{print $NF}') > tmp.txt 
    cat tmp.txt > tmp0.txt
done
cat tmp.txt | python3 plt_rel_errors.py # pipe errors into plt_rel_errors.py