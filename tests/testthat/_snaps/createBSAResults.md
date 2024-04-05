# createBSAResults works

    Code
      createBSAResults(1L, bsae, "g", window_size = 5000)
    Output
      BSAResults with 5 rows and 8 columns
                           SLB0025.vs.SLB0021 delta_alt_frequency_smoothed         g
                                    <numeric>                    <numeric> <numeric>
      CP022321.1:222-G-A           -0.0558211                   -0.0526529   1.45917
      CP022321.1:558-C-G           -0.0374720                   -0.0528294   1.00777
      CP022321.1:1443-C-T          -0.0651915                   -0.0532942   3.08363
      CP022321.1:10528-A-T         -0.0517636                   -0.0880162   0.29147
      CP022321.1:12097-T-A         -0.1361627                   -0.0959485   3.91516
                            g_smooth    filter      pvalue bh_adj_pvalue    qvalue
                           <numeric> <logical>   <numeric>     <numeric> <logical>
      CP022321.1:222-G-A     1.83249     FALSE 5.75892e-01   0.575892124        NA
      CP022321.1:558-C-G     1.85515     FALSE 4.64968e-01   0.575892124        NA
      CP022321.1:1443-C-T    1.91484      TRUE 2.09656e-01   0.524139440        NA
      CP022321.1:10528-A-T   1.84798     FALSE 5.00000e-01   0.575892124        NA
      CP022321.1:12097-T-A   2.18856      TRUE 6.05738e-05   0.000302869        NA

