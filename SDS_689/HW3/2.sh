# test different LSTM widths and plot the MSE loss vs. epoch
echo -n "" > tmp.txt # to pipe into plotting script
widths='10 20 30 40 50'
for w in ${widths}
do
    echo $w | python3 2.py > tmp0.txt # run with width w
    grep -i ms/step tmp0.txt | awk '{print $NF}' > tmp1.txt # get losses vs. epoch
    # add relative L2 error for testing
    tail -1 tmp0.txt | awk '{print $NF}' >> tmp1.txt 
	paste tmp.txt tmp1.txt > tmp0.txt # accumulate list of losses in tmp.txt
    cat tmp0.txt > tmp.txt
done
# add widths to end of tmp.txt then pipe into 2_plot.py for plotting
echo $widths >> tmp.txt 
cat tmp.txt | python3 2_plot.py # pipe tmp.txt into 2_plot.py