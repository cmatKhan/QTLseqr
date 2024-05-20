# createBSAResults works

    Code
      createBSAResults(1L, bsae, "g", window_size = 5000)
    Output
      BSAResults with 5 rows and 13 columns
        min_depth delta_snp_index delta_snp_index_smoothed     CI_05     CI_95
        <integer>       <numeric>                <numeric> <numeric> <numeric>
      1       101      -0.0558211               -0.0526529 -0.138614  0.148515
      2       188      -0.0374720               -0.0528294 -0.122340  0.117021
      3       167      -0.0651915               -0.0532942 -0.119760  0.125749
      4        37      -0.0517636               -0.0880162 -0.216216  0.216216
      5        39      -0.1361627               -0.0959485 -0.205128  0.205128
           CI_025    CI_975         g g_smoothed    filter      pvalue bh_adj_pvalue
        <numeric> <numeric> <numeric>  <numeric> <logical>   <numeric>     <numeric>
      1 -0.168317  0.168317   1.45917    1.83249     FALSE 5.75892e-01   0.575892124
      2 -0.143617  0.143617   1.00777    1.85515     FALSE 4.64968e-01   0.575892124
      3 -0.143713  0.143713   3.08363    1.91484      TRUE 2.09656e-01   0.524139440
      4 -0.243243  0.243243   0.29147    1.84798     FALSE 5.00000e-01   0.575892124
      5 -0.256410  0.230769   3.91516    2.18856      TRUE 6.05738e-05   0.000302869
           qvalue
        <logical>
      1        NA
      2        NA
      3        NA
      4        NA
      5        NA

