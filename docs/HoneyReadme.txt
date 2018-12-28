Honey Documentation

== CONTENTS == 

I.   Using This README
II.  Requirements
III. Installation
IV.  Quick Start
V.   Toy Data 
VI.  Pie Details
VII. Tails Details
VIII.Spots Details
IIX. Interpreting Results
IX.  Extras
IX.  FAQ

== I. Using This README ==

 Toy Data: 
    Provided with this distribution of Honey is a toy example 
    inside of docs/honeyExample directory. Use this once you've 
    setup Honey to test that everything is working as expected 
    and to become familiar with the software. 

 Commands: 
    All commands are presented in the format
    > commandToExecute 
    where commandToExecute would be the actual command. 

== II. Requirements ==
    
 *  samtools        0.1.7
 *  blasr           1.3.1
 *  python          2.7
 *  h5py            2.3.1
 *  pysam           0.8
 *  numpy           1.6
 *  pbdagcon        See https://github.com/PacificBiosciences/pbdagcon 
 		 	stable as of Dec 5 Commit
		 	a5f71e	709ea3590f058cfa2d2eba77fdeedf7395)


 Note: pbdagcon isn't completely required. There is a built in module that
 performs a similar function of pbdagcon, just slowly and less accurate. See
 Honey.py spots -h for details

== III. Installation ==
 
 1) Ensure all of the above requirements are in your environment.
    I recommend using python virtual environments and pip to 
    install all the packages.
 
 2) Edit setup.sh and change $SWEETPATH to the full directory where
    you've placed the package.
 
 3) To automatically place the software into your environment, add
    > source <path to>/setup.sh
    to your .bash_profile

  Be sure to source your .bash_profile (or just setup.sh) before 
  using Honey
  
  4) blasr and pbdagcon must be in your environment's $PATH

== IV. Quick Start == 
 
 Here's a quick start guide to the Honey specific steps. For more 
 details on each step in the entire process, see Section V
 below. Otherwise, check out the --help for available parameters.
 
 Beginning with a bam with mapped PacBio reads, you'll use Honey.py to
 execute the stages 'pie', 'tails', 'spots'.

 1) Honey.py pie 
    Map reads using blasr (optional) and then extract the soft-clipped tails
    and attempt to remap them.
    
 2) Honey.py tails
    Cluster the tail-mapping information to make genomic breakpoints
    
 3) Honey.py spots
    Look for genomic variants within the span of reads. 
 
 4) Honey.py force
    Look for read evidence to support a predicted structural variant
 
 5) Honey.py asm
    Assemble reads around a structural variant and map the contigs back to 
    the provided reference.
 
 6) Honey.py cpxres (beta)
    Given a set of breakpoints from tails, trace the breaks to resolve
    the full variant structure.

== V. Toy Data == 

 Feel free to use the toy data in docs/honeyExample to test your installation
 and become familiar with the Honey procedure. 
 See honeyExample/workflow.sh for the full set of commands.

== VI. Pie Details == 
 
 Command:
 > Honey.py pie inputReads.fastq reference.fasta
 
 Description:
 Pie is a wrapper around blasr that allows you to start from filtered subreads
 If you specify your input as fasta|fastq|fofn files, Pie will perform
 the primary alignments of your reads. Recommended parameters are built
 into the defualts, but you can customize them as needed.
 
 Next, pie will extract  unmapped, soft-clipped tails.
 Either the .sam/bam created in step one or if you specify your input to be
 a .sam/bam file, some number of your reads will have unmapped,
 soft-clipped tails. Pie will extract those tails, remap them, and then
 create a consolidated .sam/bam.

 Parameter Details:
 > Honey.py pie --help 
 Most of these parameters are very simple. If you use --chunk, you'll only
 create chunks of commands which are printed to stdout and can be used to
 process several inputs separately. 
 
 Follow Up:
 Once your full tails.bam is created, you should merge all your results (if
 you used --chunks) into a single sam/bam.
 
 If your merging procedure didn't already create a sorted.bam for you, use
 sam2bam to convert a single .sam to a sorted .bam and index your alignments.
 >  sam2bam reference.fasta mapping.tails.sam

== VII. Tails Details ==

 Command:
 > Honey.py tails mapping.tails.bam
 
 Description:
 This will take all of your split reads and create clusters of reads that
 split at spots within --buffer. If more than --minNumReads and --minZMWs cluster
 within any point, we'll report that cluster. The minMapq parameter will throw
 out any read with either p-i or i-e tail combination with a single alignment
 lower than the parameter. The fastq parameter will create a tar.gz where
 each read in each cluster is output into a single fastq. The verbose
 parameter will create a full report of every single read's annotation for
 every single cluster.

 
== VIII. Spots Details ==  

 Command:
 > Honey.py spots mapping.tails.bam

 Program Description:
   Spots looks for regions of intra-read discordance to discover SVs.
 The three steps in the procedure are ErrorCounting, SpotCalling, and 
 Consensus Calling. 
   First, every read is piled-up against the reference and the number of 
 insertions and deletions as well as coverage is counted and stored in 
 a hon.h5 file. In order to reduce the number of error bases that contribute 
 to the pileup, only indels greater than or equal to --minIndelErr are counted
 in the hon.h5, anything smaller is ignored. 
   Second, spots are called by looking for regions where the total number of
 insertions or deletions is greater than or equal to --threshold. Each of
 these regions are then filtered by fetching all the reads within the region
 and finding every read with >= --minIndelSize total indel bases within the
 region. Only spots that have >= --minErrReads will pass the filtering step.
   Finally, the spots are passed passed into a consenus calling algorithm so
 that de novo local assembly is created. All reads that were reported as 
 supporting the variant during the SpotCalling phase are extracted and trimmed
 down to their boundaries within the region +/- --buffer. A backbone read is
 chosen and all other supporting reads are mapped to it. The resulting
 alignment is fed into --consensus so that a high quality sequence is created.
 This sequence is the remapped to the region +/- --buffer and all SVs >=
 --minIndelSize are reported.

== VI. Output Format ==

 Honey spots:
    .spots output -- Your variant calls with the format
    * CHROM       Reference entry where variant occurs
    * START       The exact start of the variant
    * END         The exact end of the variant
    * TYPE        Variant type. One of INS (insertion) or DEL (deletion)
    * SIZE        Size of the variant
    * INFO        More information associated with the calls
    
    .h5 Output -- Contains each base-alignment type and coverage for each
    position in your reference. Use this if you wish to recall spots
    with different parameters (see recallSpots.py --help)
    
    .reads output -- If the --reads parameter is specified, a file is created
    that reports the variants in the .spots format and each line after has var
    or ref followed by a read name. This represents which reads held the
    variant (var) or which reads didn't (ref)

 Honey tails:
    .tails output -- Your variant calls with the format
    * id        Unique identifier of the call
    * chrKey    Key for what two chromosomes the breakpoint pair hits
    * uRef      First breakpoint location's reference
    * uBreak    First breakpoint location's coordinate
    * uMapq     Average mapping quality of reads that support the
                first breakpoint
    * dRef      First breakpoint location's reference
    * dBreak    First breakpoint location's coordinate
    * dMapq     Average mapping quality of reads that support the
                first breakpoint
    * remainSeq Average amount of sequence left unmapped between breakpoints
    * annot     Honey's predicted structural variant type. Note this
                isn't 100% accurate. Insertion (INS), deletion (DEL), or
                translocation (TLOC) are possible values.
    * numReads  Number of reads supporting structural variant
    * numZMWs   Number of ZMWs producing reads that support structrual
                variant
    * evidence  Semi-colon delimited list of read orientations around
                 breakpoint. See Section VII for details
    
    Honey tails can produce two extra files with data that contains
    infomation which may help validation and annotation of tail
    structural variants. The two flags to create these files are:
    1) --fastq 
       Inside of a .tgz file, write a .fastq and .bam for each cluster
       identified as clu<id>.[bam|fastq] where <id> is column 1 of
       the tgraf output. These two files contain the full sequence and
       alignment information of each initial read that supported the
       structral variant.
    2) --verboseFile
       Write a .verbose with full details of every reads' alignment
       pair that supported the structural variant
       Each cluster begins with an identifier line
       
       ##Cluster <id> - <chrKey> 
	
       followed by 1+ lines fore each read's alignment. The coloumns
       are:
       
       #uRef uBreak uMapq dRef dBreak dMapq remainSeq break annot readName
       
       These are the same as what's reported above except on a
       per-read basis and the read's name is repoted as the final
       column.
           
== VII. Interpreting Results ==
  This section hasn't been updated in a while....

  This section contains a tutorial on how to look at the data as well as
  some extra detail and advanced techniques for using Honey! To get the most
  out of this section, be sure you've run through Section V and created all of
  the results for the Toy Data
  
  = Spots Results =

  If you've followed along with the toy example, you should have found
  the following variants in the mappingFinal.hon.spots file
  
  #CHROM	OUTERSTART	START	INNERSTART	INNEREND	END	OUTEREND	TYPE	SIZE	INFO
	<removed>

  = Manual Inspection =
  
  Additionally, you can go into the bam itself to look at the alignments
  around the deletion breakpoints.
      > samtools tview mappedFinal.bam lambda_modified.fasta
     
  Press the '/' key to pull up the "Jump to Location" dialogue. Then type in
  "lambda_NEB3011:35000" then 'enter'. Move the cursor a little to the left (press 'j')
  and you should see something like this:
  
     34981          34991        35001	   35011      35021	35031	35041     35051   
  *C*ATTA*G*TG*AG*TTGA*TTG*AG*CTT*GGAATCAG*GAAGCTACGTT*CAACTCGA*****CTTATAAGGCGG*TGCCAGATG
   . .... . .. .. .... ... .. ... ..KK.KK. .KK.K.KT... KKKK.K.K     K..K.KK..K.. ..KKK.K..
  *,*,,,,*,*,,g,,*,,,,*,*,*,,*,,,*,,
  *.*....*.*..*..*....*...*..*...*..******************************************************
  *.*....*.*..*..*....T...*..*...*********************************************************
  *.*....*.**.*..*....*...*..C..**********************************************************
  ***....G.*..*..*....*...*..*...*
  *,*,,,,*,****,,*,*,,*,,,*,,*,,,*,,
  *,c,,,,*,*,,*,,g,,,,*,,,*,,*,,,t,*****,,t,**********************************************
  *,*,,,,*,g,,*,,*,,,,*,,,g,,*,,,*,*******************************************************
  *.*....*.*..*..*....*.C.*..**..*..************.GT...GA...*.*.AAAAC...********.T.T......*
  C.*....*.*..*..*....*...*..C...*********************************************************
  
  The top row is your position, second row is your reference sequence, third
  row is your consensus sequence, followed by a number of reads' alignments.
  As you can see, at about 35000, several reads either stop mapping or begin
  placing many deletion errors in the alignment. Because of the high error
  rate and alignment ambiguities, the exact point at which to begin the
  deletion varies on a per-read basis. This is straight forward in our
  simulated reference, but in real structural variants, there's usually
  microhomology around the breakpoints that ambiguity. Since our breakpoint
  here is so well definined, Honey predictes the exact breakpoint perfectly at
  position 35000. 
  
  = Recall Spots =
  
  A very time consuming part of spots calling is counting the errors. That's
  why Honey.py saves a hon.h5 file which you can provide to Honey.py spots and
  skip the time consuming step. You can change the parameters as you wish to
  see how they affect your results.
  
  > Honey.py spots mappingFinal.bam ~/reference.fasta --hon mappingFinal.hon.h5
  
  Be careful not to overwrite results that you'd like to keep.
  
  = Tails Results =
  !! WARNING -- This documentation is also incomplete
  
  Check out your tails results in mappingFinal.hon.tails and you'll
  see a total of 5 structural variants called (ids 0-3).

  ID=0 is simply evidence of lambda's circular genome. This is
  annotated as an Insertion because technically you can insert the
  sequence from the end of the reference at the beginning of the
  reference and still correctly represent your sample's genome.
  (However you wouldn't remove the INS SV call)

  ID=1 is also evidence of lambda's circular genome. However! this
  data indicates that a portion of the lambda molecules don't contain
  the first appx 900bp of the reference around their orgin. This is
  corectly called as a separate event.

  ID=2 Shows evidence of a 140bp (remainSeq) INS.
  ID=3 is an inversion. The simulated inversion was placed from
       9000-12000bp in the reference.
  ID=4 Shows a 1003bp deletion between 20000-21003, which is close to
       the exact coordinates of 20000-21003. See Manual Inspection
       above to see why the breakpoints may have been off by 3bp.
  
  = Tails Evidence =
  
  The final column in the Tails output shows how each read broke apart
  across the uBreak and dBreak. Every read starts with an initiail [i]
  alignment and any softclipped tails from the 5' of [i] that map
  create prologs [p] and from the 3' of [i] create epilogs [e]. The
  way to read a piece of evidence is to trace the read. In our sample
  genome - the sequence is always 
  
  	->p=i->i=e->
   
  Tracing down the 'pie' shows you the structure of the variant when
  represented in the reference space. 
  
  A basic example would be mapping a deletion in your sample relative to your
  reference. When taking the 'pie' structure above, we may see the orientation
 
 	->i=e->
  
  If laying this evidence out in the reference space, this would look similar
  to the following:   (ATCG) represents reference bases

	    ->i          e->
  	ACAATAGAGAAGCGACTTAGCTAGCAG
  
  Here, the sequence 'AGAAGCGACT' has been deleted in our sample.
  Another possible orientation for this deletion is ->p=i-> if our initial
  alignment hit the downstream breakpoint on the direct strand.
  
  For another example, if you look at Tail
  ID=1, you'll see
  
        e->=->i
	i<-=<-e
	i->=->p
	p<-=<-i
	
  This says that if you move on the direct strand to the dBreak (->i)
  in the reference, you'll pick back up at the uBreak and move
  downstream on the direct strand (e->). The equal sign (=) in the middle
  represents that the two pieces of the alignment are on the same
  strand. If the pieces mapped to the compliment strand, the arrows
  would be backwards (<-) and if the pieces mapped on different
  strands, the equal sign would become a percent sign (%).
  
  
  = Important Tail Parameters =
 
  The best way to learn about the tail parameters are just to mess
  with them. Rerun Honey tails
  
      > Honey.py tails mappingFinal.bam -B 1000 -f -v -o tail2
     
  You'll see that we now make 17 structural variants! Using a buffer
  of 1000bp meant that our two structural variants that showed
  evidence of the circular genome has pooled because the reads all had
  the same orientation around the breakpoints and the breakpoints were
  within the specified distance. This is an important consideration
  when trying to find variants that border large repeats. If it's
  possible for the read to interrupt it's mapping over a large area, a
  larger buffer will collect all of the reads that support the event.
  However! You may incorrectly pool the breakpoints of two separate
  structural variants that happen to be in the same region within the
  buffer's size distance.
  
  = Missed Adapters =
  
  For the other events created with the modified parameters, these are
  evidence of a missed adapter. The way to tell is if the two
  breakpoints are near each other (usually within ~50bp or ~100bp) and
  the sample sequence switches strands and begins mapping back over
  itself (->%<- or <-%->). One way to possibly filter these out is to
  increase the number of reads and number of ZMWs required to support
  an event before it's reported. The below command will filter all
  but one of our missed adapter events.
  
      > Honey.py tails mappingFinal.bam -B 1000 -f -v -o tail3 -b 6 -z 6
  
  These can be identified as missed adapters by then looking at the
  verbose output provied in tail3.verbose. As you can see for the
  reads that support tail ID=3, they all look independent missed
  adapters. The problem is that they all occured close enough that
  they fell with the buffer. 

== VIII. Extras ==
  Change your calls into a .bed by using spotToBed.py and
  tailToBed.py

== IX. FAQ ==
  
  * Who can I report bugs to or ask questions?
    - Please report your issues to the sourceforge ticketing system.

  * I get the pyvcf error 
  	"error: can't copy 'vcf/cparse.c': doesn't exist or not a regular file"
    - To fix this, manually compile vcf/cparse.pyx using cython
      	cython cparse.pyx
	
