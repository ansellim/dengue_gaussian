Initialize a git repository, and Git add and commit the relevant files, then make the following changes:
  1. Remove wolbachia coverage and NPI stringency as covariates. Only keep temperature and rain.                 
  2. Since i am modeling cases rather than infections, i should model for a delay in detection of cases. What  
  would be appropriate based on the literature? Either i include this term in the convolution, or i should       
  model infections ...Also do mention somewhere that i am assuming 100% Case ascertainment rate ...              
  3. Try reducing the amount of noise in the Gaussian process and see whether it will affect the variance        
  decomposition.                                                                                                 
  4. The focus is on developing an early warning system for increases in Rt using serotype switching (changes    
  in serotype distribution, including entropy rise and fall). Are increases in Rt contemporaneous with,          
  preceded by, or followed by, serotype switching? Also examine the hypothesis that in the weeks prior to        
  serotype switching, we can observe an increase in Rt, and how soon we can observe this rise. Maybe plot the    
  entropy rise and fall around each serotype switching changepoint to see it more clearly.