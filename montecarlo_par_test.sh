TEST_RESULTS_FILE="montecarlo_par_test_results.txt"
TIME_RESULTS_FILE="montecarlo_par_time_results.txt"

rm "$TEST_RESULTS_FILE" 2>/dev/null
rm "$TIME_RESULTS_FILE" 2>/dev/null

for ((N=1; N<=1000000000; N*=10))
do
    OUTPUT=$(./gauss2d_montecarlo_integral_par "$N")
    
    echo -n "$N Hiba: " 
    echo "$OUTPUT" | grep Hiba | awk '{print $2}' | tee -a "$TEST_RESULTS_FILE"
    
    echo -n "Idő: "
    echo "$OUTPUT" | grep Idő | awk '{print $2}' | tee -a "$TIME_RESULTS_FILE"

    echo ""
done