#!/usr/bin/env python3
##################################################################
## Jose F. Sanchez                                              ##
## Copyright (C) 2019 - 2023 Lauro Sumoy Lab, IGTP, Spain       ##
##################################################################
'''
Generates sample quality control at raw reads, assembly or annotation level.
'''
## useful imports
import time
import os
import concurrent.futures
import pandas as pd
from termcolor import colored

## import my modules, scripts, config
from BacterialTyper.scripts import fastqc_caller
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.modules import help_info
from BacterialTyper import __version__ as pipeline_version

from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info

################################################
def run_QC(options):

    """
    Main function of the QC module.
    
    It creates a quality check for each sample:
        - it can be generated at the raw reads or trimmed reads, using fastqc and summarizes using MultiQC
        - it can be generated using assembled contigs or annotated proteins to check the assembly completeness
    
    .. seealso:: This function depends on other BacterialTyper and HCGB functions called:
    
    - :func:`BacterialTyper.scripts.BUSCO_caller.print_help_BUSCO`
    
    - :func:`BacterialTyper.scripts.multiQC_report.multiqc_help`
    
    - :func:`HCGB.sampleParser`
    
    - :func:`HCGB.functions.aesthetics_functions`
    
    - :func:`HCGB.functions.time_functions`
    
    - :func:`HCGB.functions.main_functions`
    
    - :func:`HCGB.functions.file_functions`
    
    .. include:: ../../links.inc         
    
    """
    
    ## init time
    start_time_total = time.time()

    #################################
    ## show help messages if desired    
    #################################
     
    # if any help_flag provided will print and exit
    if (help_info.help_info(options) == 1):
        raise SystemExit() 
        
    ## debugging messages
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False    
    
    ### set as default paired_end mode
    if (options.single_end):
        options.pair = False
    else:
        options.pair = True
    
    ## set main header
    HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
    HCGB_aes.boxymcboxface("Quality check")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## Project mode as default
    if (options.detached):
        options.project = False
        outdir = os.path.abspath(options.output_folder)
    else:
        options.project = True
        outdir = input_dir        
    
    ### option
    if (options.raw_reads):
        ## get files
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, 
                                                      "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], 
                                                      options.debug)
        fastqc(pd_samples_retrieved, outdir, options, start_time_total, "raw", Debug)
        submodule_name = "qc_raw_reads"
    elif (options.trim_reads):
        ## get files
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
        fastqc(pd_samples_retrieved, outdir, options, start_time_total, "trimmed", Debug)
        submodule_name = "qc_trimm_reads"
    elif (options.assembly):
        pd_samples_retrieved = BUSCO_check(input_dir, outdir, options, start_time_total, "genome")
        submodule_name = "qc_assembly"
    elif (options.annotation):
        pd_samples_retrieved = BUSCO_check(input_dir, outdir, options, start_time_total, "proteins")
        submodule_name = "qc_annot"
        
    print ("\n*************** Finish *******************")
    HCGB_time.timestamp(start_time_total)
    
    ## samples information dictionary
    samples_info = {}
    samples_frame = pd_samples_retrieved.groupby('new_name')
    for name_tuple, grouped in samples_frame:
        name = name_tuple[0]
        samples_info[name] = grouped['sample'].to_list()
    
    ## dump information and parameters
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"qc", "time":time.time(),
                "BacterialTyper version":pipeline_version,
                'sample_info': samples_info }
    
    HCGB_info.dump_info_run(info_dir, submodule_name, options, runInfo, options.debug)
    
    ## dump conda details
    HCGB_info.dump_info_conda(info_dir, submodule_name, options.debug)
    
    print ("\n+ Exiting QC module.")
    return()

################################################
def fastqc(pd_samples_retrieved, outdir, options, start_time_total, name_analysis, Debug):
    """
    Function of the QC module to produce results for FASTQC.
    
    It creates a FASTQC quality check for each sample: it can be generated at the raw reads 
    or trimmed reads, It uses fastqc and summarizes using MultiQC.
    
    This function is called from the run_QC function in this QC module and from
    the run_trim() function from the trim module after trimming reads.
    
    It controls the read file inputs provided and calls fastqc_caller.run_module_fastqc()
    Then it generates a report using MultiQC.
    
    It uses variables name_analysis for creating different analysis: raw, trimmed, etc
    
    :param pd_samples_retrieved: Dataframe containing samples information generated using HCGB.sampleParser
    :param outdir: Absolute path file to store results
    :param options: DIctionary of options as retrieved from BacterialTyper
    :param start_time_total: Timestamp for starting time
    :param name_analysis: Type of analysis
    :param Debug: Print debug messages or not
    
    :type pd_samples_retrieved: pd.dataFrame
    :type outdir: string
    :type options:  dictionary
    :type start_time_total: timestamp
    :type name_analysis: string
    :type Debug: boolean
    
    """
    
    HCGB_aes.boxymcboxface("FASTQC Quality check for samples")
    
    ## debug message
    if (Debug):
        print (colored("\n**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)
        print ("\n")

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    
    ## if not project, outdir contains the dir to put output
    ## in this case, in some other cases might not occur    
    if not options.project:
        HCGB_files.create_folder(outdir)
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "fastqc_" + name_analysis, options.debug)
    
    print ("+ Checking quality for each sample retrieved...")
    start_time_partial = start_time_total
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["name"])

    ## optimize threads
    name_list = set(pd_samples_retrieved["name"].tolist())
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        HCGB_aes.debug_message("options.threads: " + str(options.threads), "yellow")
        HCGB_aes.debug_message("max_workers: " + str(max_workers_int), "yellow")
        HCGB_aes.debug_message("threads_job: " + str(threads_job), "yellow")

    ## send for each sample
    print ("+ Calling fastqc for samples...")    
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
        commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, outdir_dict[name[0]], sorted( cluster["sample"].tolist() ), name[0], threads_job): name[0] for name, cluster in sample_frame }
        
        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("+ FASTQC for samples has finished...")    
    
    ## functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)

    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = HCGB_files.create_subfolder("report", outdir)

        ## get subdirs generated and call multiQC report module
        givenList = []
        print ("+ Detail information for each sample could be identified in separate folders:")
        
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)

        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")
        
        fastqc_report = HCGB_files.create_subfolder("FASTQC", outdir_report)
        fastqc_final_report = HCGB_files.create_subfolder(name_analysis, fastqc_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "FASTQC", fastqc_final_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %fastqc_final_report)

    return()


################################################
def BUSCO_check(input_dir, outdir, options, start_time_total, mode, Debug=False):
    
    """
    Function of the QC module to produce results for BUSCO.
    
    It creates a BUSCO quality check for each sample: it can be generated using 
    assembled contigs or annotated proteins to check the assembly completeness

    This function is called from the run_QC function in this QC module and from
    the run_assemble() & run_annot() function from the assemble and annot modules.
    
    It controls the fasta file inputs provided and calls BUSCO_caller.BUSCO_call()
    Then it generates a report using MultiQC.
    
    It uses variables mode for creating different analysis: genome or proteins
    
    :param input_dir: Folder containing files to check: assembled genomes or protein fasta sequences.
    :param outdir: Absolute path file to store results
    :param options: DIctionary of options as retrieved from BacterialTyper
    :param start_time_total: Timestamp for starting time
    :param mode: Type of analysis
    :param Debug: Print debug messages or not
    
    :type input_dir: pd.dataFrame
    :type outdir: string
    :type options:  dictionary
    :type start_time_total: timestamp
    :type mode: string
    :type Debug: boolean
    
    """

    HCGB_aes.boxymcboxface("BUSCO Analysis Quality check")

    ## absolute path for in & out
    database_folder = os.path.abspath(options.database)

    ## get files and get dir for each sample according to mode
    if mode == 'genome':
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, 
                                                      "assembly", ["fna"], 
                                                      options.debug)

        if not options.project:
            outdir = HCGB_files.create_subfolder("assembly_qc", outdir)

        if options.debug:
            print ("** DEBUG: pd_samples_retrieved")
            print (pd_samples_retrieved)
        
        BUSCO_outdir_dict = HCGB_files.outdir_project(outdir, options.project, 
                                                pd_samples_retrieved, "assemble_qc", 
                                                options.debug)

    elif mode == 'proteins':
        pd_samples_retrieved = sampleParser.files.get_files(options, outdir, 
                                                      "annot", ["faa"], 
                                                      options.debug) ##

        if not options.project:
            outdir = HCGB_files.create_subfolder("annot_qc", outdir)

        if options.debug:
            print ("** DEBUG: pd_samples_retrieved")
            print (pd_samples_retrieved)
            
        BUSCO_outdir_dict = HCGB_files.outdir_project(outdir, options.project, 
                                                pd_samples_retrieved, "annot_qc", 
                                                options.debug)

    ## add column to dataframe
    pd_samples_retrieved['busco_folder'] = ""
    for index, row in pd_samples_retrieved.iterrows():
        pd_samples_retrieved.at[index, 'busco_folder'] = BUSCO_outdir_dict[ row['name'] ] 

    ## debug message
    if (options.debug):
        HCGB_aes.debug_message("df_samples_busco", 'yellow')
        print (pd_samples_retrieved)
        
        HCGB_aes.debug_message("BUSCO_outdir_dict", 'yellow')
        print (BUSCO_outdir_dict)

    ## Check each using BUSCO
    database_folder = os.path.abspath(options.database)
    BUSCO_Database = HCGB_files.create_subfolder('BUSCO', database_folder)
    if not os.path.exists(BUSCO_Database):
        HCGB_files.create_folder(BUSCO_Database)

    ## call
    (dataFrame_results, stats_results) = BUSCO_caller.BUSCO_call(options.BUSCO_dbs, 
                                                              pd_samples_retrieved, 
                                                              BUSCO_Database, options.threads, mode, Debug)
    
    ## debug message
    if (options.debug):
        HCGB_aes.debug_message("dataFrame_results", 'yellow')
        HCGB_main.print_all_pandaDF(dataFrame_results)
    
        HCGB_aes.debug_message("stats_results", 'yellow')
        HCGB_main.print_all_pandaDF(stats_results)
        
    
    ## functions.timestamp
    print ("+ Quality control of all samples finished: ")
    HCGB_time.timestamp(start_time_total)
    
    ## multiqc report plot
    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report BUSCO plot.")
        outdir_report = HCGB_files.create_subfolder("report", outdir)

        ## get subdirs generated and call multiQC report module
        print ("+ Detail information for each sample could be identified in separate folders.")
        
        ## name folder according to mode
        if mode == 'genome':
            BUSCO_report = HCGB_files.create_subfolder("BUSCO_assembly", outdir_report)
        elif mode == 'proteins':
            BUSCO_report = HCGB_files.create_subfolder("BUSCO_annot", outdir_report)

        ## generate plots
        print ("+ Generate summarizing plots...")
        BUSCO_caller.BUSCO_plots(dataFrame_results, BUSCO_report, options.threads, Debug)    
        print ('\n+ Check quality plots in folder: %s' %BUSCO_report)

        ## TODO: Fix this chunk of code
        ##    Parse BUSCO statistics in dataframe (stats_results) for discarding samples if necessary
        ##    given a cutoff, discard or advise to discard some samples

        ### print statistics
        stats_results.to_csv( os.path.join(BUSCO_report, "BUSCO_stats.csv"))
        name_excel = os.path.join(BUSCO_report, "BUSCO_stats.xlsx")
        
        ## save in excel
        with pd.ExcelWriter(name_excel, engine="xlsxwriter", engine_kwargs={"options": {"nan_inf_to_errors": True}}) as writer:
            stats_results.to_excel(writer, sheet_name="BUSCO statistics")    
        
        print ('\n+ Check quality statistics in folder: %s' %BUSCO_report)
    
    return(dataFrame_results, pd_samples_retrieved)
