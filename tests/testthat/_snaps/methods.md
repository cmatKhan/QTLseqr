# depth getters all samples works

    Code
      totalDepth(bsae)
    Output
                           SLB0021 SLB0025
      CP022321.1:222-G-A       290     101
      CP022321.1:558-C-G       313     188
      CP022321.1:1443-C-T      292     167
      CP022321.1:10528-A-T      59      37
      CP022321.1:12097-T-A      39      58

---

    Code
      altDepth(bsae)
    Output
                           SLB0021 SLB0025
      CP022321.1:222-G-A        65      17
      CP022321.1:558-C-G        70      35
      CP022321.1:1443-C-T       61      24
      CP022321.1:10528-A-T      19      10
      CP022321.1:12097-T-A       8       4

# depth getters sample select works

    Code
      totalDepth(bsae, sample_set = "SLB0021")
    Output
                           SLB0021
      CP022321.1:222-G-A       290
      CP022321.1:558-C-G       313
      CP022321.1:1443-C-T      292
      CP022321.1:10528-A-T      59
      CP022321.1:12097-T-A      39

---

    Code
      altDepth(bsae, sample_set = "SLB0025")
    Output
                           SLB0025
      CP022321.1:222-G-A        17
      CP022321.1:558-C-G        35
      CP022321.1:1443-C-T       24
      CP022321.1:10528-A-T      10
      CP022321.1:12097-T-A       4

# snpIndex works

    Code
      actual
    Output
                             SLB0021
      CP022321.1:222-G-A   0.2241379
      CP022321.1:558-C-G   0.2236422
      CP022321.1:1443-C-T  0.2089041
      CP022321.1:10528-A-T 0.3220339
      CP022321.1:12097-T-A 0.2051282

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

# aggregate works

    Code
      actual
    Output
      <simpleMessage in message("You will may want to re-do the `comparisons` frame and add ",     "variables to the colData"): You will may want to re-do the `comparisons` frame and add variables to the colData
      >

