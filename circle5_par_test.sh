TEST_RESULTS_FILE="circle5_par_test_results.txt"
TIME_RESULTS_FILE="circle5_par_time_results.txt"

rm "$TEST_RESULTS_FILE" 2>/dev/null
rm "$TIME_RESULTS_FILE" 2>/dev/null

for ((N=1; N<=10000; N*=10))
do
    N_PHI=$((N * 6))
    
    OUTPUT=$(./gauss2d_circle_integral_par 5 "$N" "$N_PHI")
    
    echo -n "$N Hiba: " 
    echo "$OUTPUT" | grep Hiba | awk '{print $2}' | tee -a "$TEST_RESULTS_FILE"
    
    echo -n "Idő: "
    echo "$OUTPUT" | grep Idő | awk '{print $2}' | tee -a "$TIME_RESULTS_FILE"

    echo ""
done