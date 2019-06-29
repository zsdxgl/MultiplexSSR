Manual of multiplex PCR

1. About multiplex PCR

    The aim of the pipeline is to develop multiplex SSR-PCRs with resequencing data. This pipeline includes two scripts, MultiplexSSR.pl and random.pl. MultiplexSSR.pl takes resequencing data as input to develop multiplex SSR-PCRs. random.pl takes SSRs in vcf format as input to access the saturation of allele number, range, maximum position and minimum position with the increase of individual number for each SSR.

2. MultiplexSSR.pl packages six programs, including: 

   1) Tandem Repeats Finder		http://tandem.bu.edu/trf/trf.html
   2) Lobstr				      		http://lobstr.teamerlich.org/
   3) BWA				          		http://bio-bwa.sourceforge.net/
   4) SAMtools		      			http://samtools.sourceforge.net/
   5) ePCR					        	https://launchpad.net/ubuntu/+source/epcr/
   6) Multiplx			      		http://bioinfo.ut.ee/?page_id=167


3. MultiplexSSR.pl includes six steps:

   1) The tandem repeats are first detected using Tandem Repeats Finder. 
   2) The SNPs and Indels are called with BWA-MEM and SAMtools, and SSRs are called with Lobstr.
   3) The flanking sequences of 500 bp in each side of SSRs are extracted. 
   4) The primers are designed with Primer3 with SNPs, Indels and tandem repeats in reference are masked.
   5) The primers are filtered out with ePCR that match more than two locations in reference. 
   6) The primer pairs are grouped with Multiplx and subgrouped according to the position and range.
 
    The designed primers are subgrouped for dye label and could be directly used.


4. random.pl is dependent on three R pacakges

   1)	psych		
   2)	airquality
   3)	ggplot2


5. random.pl includes two steps:

   1)	The alleles for each SSR are counted for number, range, maximum size and minimum size. The maximum individual number is automatically detected and the minimum individual number is set five.
   2)	The relationship between number, range, maximum size and minimum size with individual number is tested and showed as graphs.


6. System requirement

    The pipeline was tested on Linux version 4.13.0-36-generic (buildd@lgw01-amd64-033) (gcc version 5.4.0 20160609 (Ubuntu 5.4.0-        6ubuntu1~16.04.9))


7. Data requirement.

   1) reference sequence in format of fasta
   2) resequencing data in the format “fq/fastq” without compression.


8. Create a list file: input.txt

        less input.txt
        test1.1.fq      test1.2.fq
        test2.1.fq      test2.2.fq
    
    
9. Run the pipeline

    The script MultiplexSSR.pl can be run in two models, the beginning-to-end model and skip model. 
    
    9.1 Run the pipeline in beginning-to-end model. Under this model, the script is run from the first step to the sixth step.
       
       MultiplexSSR.pl –t 20
    
    9.2 Run the pipeline in skip model. This model skip the first two steps with two advantages. The step of genotyping is skipped, which is not necessary for parameter adjustment and time wasting. The other sourced genotypes of SSRs and SNPs could be conveniently used. Only ePCR is dependent.
    
       MultiplexSSR.pl –s T –d 5
    
10. The results

	1) The primers for each SSR
	2) The SSR groups
    
11. Check the saturation with individual number
    
        random.pl –i SSR.vcf
