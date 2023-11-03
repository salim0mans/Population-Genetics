# Masking  

## Introduction     

Masking is an important method to artificially check the accuracy of a model. It includes intentional removal (or hiding) of an actual observed part of the data, and letting the model make prediction for the hidden parts. Eventually, We compare the observed "masked" parts with the prediction by machine, and decide the accuracy after getting parameters of a confusion matrix (TP,FP,TN,FN).  

There are different types of masking in genotype files: masking a whole chromosome, or a whole continuous region within a chromosome; multiple random regions are also possible. However, prediction depends on LD scores, which show a higher association between closer variants, and create separate blocks of LD that are handled as kind of less dependent units. Eventually, masking a whole region or a chromosome would wipe out whole LD blocks and add huge difficulty to the prediction procedure (hence, reduce its accuracy).  
One of the best solutions is random variant masking, which is removing a specified proportion of random variants (rows), which is topic of this protocol.   

One important point to consider is that even using a random exclusive sampling doesn't make the output representative of the population (the original data). Therefore, a cross-validation method could do a great deal negating any bias of that kind.     
Cross-validation considers basically repeating the masking procedure multiple times using different seeds, so at best scenario, a sufficent combinations of masked variants being tested, and contribute to a mean value that represents the whole procedure better.    

## Procedures genotype masking  

### Use shuf command    

We can take all IDs from the vcf file by bcftools query, and execute shuffle command on them.   
Then, calculate the number of variants to mask according to the masking percentage, and then apply the number of shuffled sample of IDs equal to that number and save this random sample of IDs into a new text file.     
Finally, we exclude the random sampled ids from the original vcf file using bcftools view -e     
  
*The masking code is available as bash scripts under the names of "single_file_shuf_masking.sh" and "parallel_file_shuf_masking.sh" in the "Scripts" directory.     
   
* Single file masking: requires the input of filename prefix and percentage of random masking   
* Parallel file masking: requires the use of "parallel" package putting filename list as the first argument and percentages as the second.  

$ parallel ./parallel_shuf_masking.sh ::: file1 file2 file3 ::: $(seq 5 5 25) ::: $(seq 1 10) ## 3 files masked by percentages 5 10 15 20 25, 10 times (10-fold validation)   
    
As mentioned before, it is highly recommended to perform the task multiple times 
  
----------------Now We have a bunch of masked files of different masking percentages, repeated multiple times----------------
