# createSampleComparisonFrame works

    Code
      bsae@comparisons
    Output
      DataFrame with 1 row and 2 columns
        population_1 population_2
         <character>  <character>
      1      SLB0021      SLB0025

# populationDepths works

    Code
      actual
    Output
      $population_1
      $population_1$alt
                           SLB0021
      CP022321.1:222-G-A        65
      CP022321.1:558-C-G        70
      CP022321.1:1443-C-T       61
      CP022321.1:10528-A-T      19
      CP022321.1:12097-T-A       8
      
      $population_1$ref
                           SLB0021
      CP022321.1:222-G-A       225
      CP022321.1:558-C-G       243
      CP022321.1:1443-C-T      231
      CP022321.1:10528-A-T      40
      CP022321.1:12097-T-A      31
      
      
      $population_2
      $population_2$alt
                           SLB0025
      CP022321.1:222-G-A        17
      CP022321.1:558-C-G        35
      CP022321.1:1443-C-T       24
      CP022321.1:10528-A-T      10
      CP022321.1:12097-T-A       4
      
      $population_2$ref
                           SLB0025
      CP022321.1:222-G-A        84
      CP022321.1:558-C-G       153
      CP022321.1:1443-C-T      143
      CP022321.1:10528-A-T      27
      CP022321.1:12097-T-A      54
      
      
      attr(,"class")
      [1] "PopulationDepthList"

