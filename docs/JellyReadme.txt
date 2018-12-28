Jelly Documentation -- v14.1

== CONTENTS ==

I.   Using This README
II.  Requirements
III. Installation
IV.  Quick Start
V.   Running Jelly
VI.  Extras

== I. Using This README ==

 Toy Data: 
    Provided with this distribution of Jelly is a toy example 
    inside of docs/jellyExample directory. Use this once you've 
    setup Jelly to test that everything is working as expected 
    and to become familiar with the software. 

 Commands: 
    All commands are presented in the format
    > commandToExecute 
    where commandToExecute would be the actual command. 

== II. Requirements ==

 * Blasr (https://github.com/PacificBiosciences/blasr)
    Version 1.3.1 is fully vetted as compatible with Jelly. 
    Other versions may run into problems. Use 
    > blasr -version
    to figure out what you have. Blasr must be in your environment
    path. Also, you'll need sawriter (included with blasr) installed

 * Python 2.7
    Python must be in your environment path and executable with 
    the commands:
    > python
    > /usr/bin/env python
  
 * Networkx v1.7
    NOTE!! Previous versions of PBSuite required v1.1 and anything past this
    would cause bugs. Now we are compatable with v1.7 and most likely any
    version past v1.7 up to v1.9. I have updated the code try working with
    whatever version of networkx you have installed. If there is a problem,
    the logs will hold a warning and you'll need to change your version.
    If you see anything else that's weird, file a bug, please.
    
== III. Installation ==

 1) Edit setup.sh and change $SWEETPATH to the full directory where 
    you've placed the package.

 2) To automatically place the package into your environment, add
    > source <path to>/setup.sh
    to your .bash_profile.

  Be sure to source your .bash_profile (or just setup.sh) before 
  using Jelly.


== IV. Quick Start ==

  For more details on each step in the pipeline, see Section V
  below. If, however, you'd like to just run the program do the
  following.

  1) Create your Protocol.xml 
     To run the lambdaExample dataset provided, edit the paths in 
     the <reference> , <outputDir> and the baseDir attribute in 
     the inputs tag to the full path in which lambdaExample is 
     sitting. 
     See Section V.1 for details about creating Protocol.xml
 
  2) Run each stage
     Sequentially execute each state. One stage must finish executing
     before continuing to the next. To run a stage, use the command
     > Jelly.py <stage> yourProtocol.xml 

     The stages, in order, and their descriptions are
       1. setup       Tag sequence names, find gaps, and index the reference
       2. mapping     Use blasr to map the sequences to the reference
       3. support     Identify which reads support which gaps
       4. extraction  For each gap, consolidate all reads supporting it into a local-assembly folder.
       5. assembly    Build the consensus gap-filling sequence
       6. output      Stitch the reference sequences and gap-fillling sequences together.

     To get help with Jelly, simply run Jelly.py --help.
     Or, for help with any stage, simply run 
     > Jelly.py <stage> --help 
     Descriptions of each step are in Section V.
  
  3) Passing Parameters through Jelly.py
     If you would like to pass a parameter to the stage you are running, use
     "-x". For example, when running the support stage, if you only wanted
     Jelly to attempt to fill captured-gaps (i.e. no inter-scaffold gaps), and
     you wanted to require that a read must have a minimum mapping QV of >=
     250 to support a gap, you'd use the command:
     > Jelly.py support Protocol.xml -x "--capturedOnly --minMapqv=250"
     All parameters you add need to be enclosed in double quotes after the -x

== V. Running Jelly ==

 = Pre-Processing =

    1) Initial Stats
       If you would like to get some size information about your 
       reference and/or input reads run 
       >  summarizeAssembly.py <reference.fasta> 
       >  readSummary.py <Protocol.xm 
       see summarizeAssembly.py --help and readSummary.py for 
       detail 
   
    2) Prepare your input data
       Gather the paths to your reference and all of your input 
       sequence files. If using PacBio reads, use filtered subreads
       where SMRTBell adapters have been removed.
 
       If you have a small number of very large sequence
       files and you want to speed up processing, split those 
       into several smaller files. Jelly will submit one 
       mapping/support job per sequence file.

       Every sequence file you use has the following requirements:
       * Sequences can be in a .fastq file or in a .fasta file 
         (no .fa, .fsta, etc extension.). 
       * References must be .fasta
       * Fasta files (reference or reads) must have an associated 
         qual file with the same name sitting in the same directory 
	 beside the .fasta file. 
       * Qual files should contain the Phred Scores of bases (0-93) and 
         should not be encoded (i.e. no Sanger/Solexa, only the 
         number for the score) If you do not have qualities for your 
         sequences, use
         > fakeQuals.py --help 
       * Each set of input reads need to have a unique file
         name. (e.g. filtered_subreads1.fastq, filtered_subreads2.fastq)
       * Sequence names CANNOT have spaces.
         
 = Create Your Protocol =
  
    This is by far the longest and most involved step. Once you 
    get past this, Jelly makes the rest of the workflow super 
    easy. See TemplateProtocol.xml for an idea of what a protocol 
    should look like. You can name your protocol whatever you'd like. 
    
    Below are the elements needed for a Protocol.
    * <reference> : The full path to your reference.fasta
        All of the files Jelly creates regarding your 
        reference will be placed in the same directory beside 
        the reference.fasta 

    * <outputDir> : The output tag contains the full path to 
        where Jelly will put the intermediate data for each stage 
        in the process. Placing your Protocol.xml into the outputDir 
        directory makes for excellent bookkeeping. Just copy the 
        provided TemplateProtocol.xml into the outputDir directory 
        for each project you maintain 

    * <cluster> [optional]: If you do not wish to submit jobs to 
        a cluster, do not include the <cluster> element in your 
        protocol. Otherwise, the cluster tag contains information 
        that Jelly will use to submit jobs to your cluster. The 
        cluster tag contains two elements.
    
    * <command> : This is a template that holds the structure of 
        how one submits jobs to his/her cluster. The example 
        provided in TemplateProtocol.xml is used on a MOAB job 
        management system. The command templatees has 4 REQUIRED 
        elements: 

          * CMD     - The command one uses to execute on the cluster. 
          * JOBNAME - The name to assign to the job being submitted. 
          * STDOUT  - The file that standard out will be directed to. 
          * STDERR  - The file that standard error will be directed to.

        If you have a single, large machine and not a cluster,
        can use the following command to submit all of your jobs
        he background and parrallize operations. Just be careful
        the number of jobs/resources you execute or you can 
        freeze your system.
        <command>${CMD} ${JOBNAME} 2> ${STDERR} 1> ${STDOUT} &amp;</command>

    * <nJobs> : This is the maximum number of jobs Jelly can 
        submit to the cluster. Jelly tries to submit as many jobs 
        as possible. So, if you do not specify nJobs, Jelly will 
        create a job for every single input file for any given stage. 
        For example, if you have 50 files containing read sequences 
        for mapping, and you specify nJobs = 5, Jelly will submit 
        5 jobs, each will map 10 input files. If you do not specify 
        nJobs (i.e. don't include the nJobs tag or set nJobs to 0) 
        Jelly will submit 50 jobs for mapping. 

    * <blasr> : Place all of the mapping parameter to be sent to 
        blasr here. You can be comfortable with the default blasr 
        parameters that are provied in the TemplateProtocol.xml. 
        However, if you would like to customize the parameters, be 
        sure to read the blasr documentation (blasr --help).
        
        Always specify -noSplitSubreads in your blasr parameters. 

    * <input> : Input contains information about where your input 
        data is located. First, there is the optional 'baseDir' 
        attribute. If all of your data has a common root path, 
        specifying baseDir=/That/Path will prevent redundancy in 
        the inputs. Inside of the <input> tag is the <job> tag. 
        This contains the path to each input file Jelly will use. 
        
 = Setup your files =

    >  Jelly.py setup Protocol.xml

    All of the files Jelly creates regarding your reference will
    be placed in the same directory beside the reference

 = Mapping your data =

    > Jelly.py mapping Protocol.xml 

    Remember to wait until given stage is finished before running
    the next stage. 
    Standard Error and Standard Output Logs for each step are placed next 
    to the data Jelly creates. For the mapping step, this is
    found in <outputDir>/mapping/

 = Support The Gaps =

    > Jelly.py support Protocol.xml 
   
 = Extract Useful Reads =

    > Jelly.py extraction Protocol.xml 
  
 = Assemble The Gaps =

    > Jelly.py assembly Protocol.xml 
    
    If you have access to more than one core per gap  to be 
    assembled, be sure to tell Jelly to pass the nproc parameter 
    to the assembly stage via:
    >  Jelly.py assembly Protocol.xml -x "--nproc=4"
    Where 4 can be replaced by the number of cores you're using.
    
 = Output Your Results =

    > Jelly.py output Protocol.xml 
    
    At the head of your log file, you can find information 
    about how many gaps were addressed, filled, etc. The output 
    stage collects all of your results into 3 files:
    
    <outputDir>/pbjelly.out.fasta
    <outputDir>/pbjelly.out.qual
    <outputDir>/liftOverTable.json 
    
== VI. Extras ==

  * blasrToBed.py 
    This script will convert blasr's .m4 or .m5 format into a 
    BED Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
    If you would like to visualize the alignments, I 
    reccommend using IGB ( http://bioviz.org/igb/index.html ).

  * bedToCoverageWig.py 
    Turn a bed file with alignments into a depth of coverage 
    WIG Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format6 ).
 
== VI. FAQ ==

  * Who can I report bugs to or ask questions?
    Please report your issues to the sourceforge ticketing system.

