# runs 1.py and saves outputs in 1.txt (for viewing later)
# Then, processes the output and pipes it into 1_plot.py for plotting
# the cross-entropy loss and training error as a function of epoch

cd /c/Users/jonas/Documents/SciML_class/HW3
python3 1.py > 1.txt
# extract cross-entropy loss as a function of epoch
grep -i Cross-entropy 1.txt | awk '{print $NF}' > tmp.txt 
# pipe into 1_plot.py contains three columns:
# 1) cross-entropy loss as a function of epoch
# 2) training error as a function of epoch
paste tmp.txt <(grep -i Training 1.txt | awk '{print $NF}') | python3 1_plot.py