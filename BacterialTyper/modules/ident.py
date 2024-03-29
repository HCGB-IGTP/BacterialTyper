#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Generates sample identification using KMA software and MLSTar. 
Looks for similar entries on GenBank and retrieves them.  
"""
## useful imports
import time
import os
import concurrent.futures
from termcolor import colored
import pandas as pd
import pprint

## import my modules
from BacterialTyper.scripts import KMA_caller
from BacterialTyper.scripts import kraken2_caller
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import MLST_caller
from BacterialTyper.scripts import edirect_caller
from BacterialTyper.modules import help_info
from BacterialTyper.config import set_config
from BacterialTyper.scripts import database_user
from BacterialTyper import __version__ as pipeline_version

from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info


####################################
def run_ident(options):
    """
    Main function acting as an entry point to the module *ident*.
    
    This function generates an species typification for each sample retrieved using some of the following options
    1. Kmer alignment (KMA) or kraken2 software.
    2. MLST profiles based on species identification or user provided input.
    
    Additionally, it can also search against pre-defined databases by KMA or user-defined databases.

    Arguments:
    
    
    .. seealso:: Additional information to PubMLST available datasets.
    
        - :doc:`PubMLST datasets<../../../data/PubMLST_datasets>`
    
    
    """


    ##################################
    ### show help messages if desired    
    ##################################
   
    ## if any help_flag provided will print and exit
    if (help_info.help_info(options) == 1):
        raise SystemExit() 
        
    ## init time
    start_time_total = time.time()

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

    ### KMA_caller -> most similar taxa
    HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
    HCGB_aes.boxymcboxface("Species identification")

    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## Project mode as default
    global Project 
    
    if (options.detached):
        options.project = False
        project_mode=False
        outdir = os.path.abspath(options.output_folder)
        Project=False
    else:
        options.project = True
        outdir = input_dir
        Project=True

    ## print options
    if (Debug):
        HCGB_aes.print_argparse_dict(options)

    ## get files: trimmed, assembly, annot
    pd_samples_retrieved = database_user.get_userData_files(options, input_dir)
    
    ## debug message
    if (Debug):
        HCGB_aes.debug_message("pd_samples_retrieve")
        print (pd_samples_retrieved)
    
    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        HCGB_files.create_folder(outdir)
        
    ## for each sample
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, 
                                            pd_samples_retrieved, "ident", 
                                            options.debug,  groupby_col="name_sample")    
    
    ## let's start the process
    print ("+ Generate an species typification for each sample retrieved using some of the following options:")
    print ("(1) Kmer alignment (KMA) or kraken2 software.")    
    print ("(2) Pre-defined databases by KMA or user-defined databases.")
    print ("(3) MLST profiles based on species identification or user provided input.")
    
    ## default kraken database
    if options.kraken2:
        print("Kraken2 databases:")
        if options.kraken_dbs:
            print("+ Using database provided:")
            print(options.kraken_dbs)
        else:
            print("+ Using default database: PlusPFP-8 Gb")
            options.kraken_dbs = 'pluspfp_8'

    
    ## get databases to check
    retrieve_databases = get_options_db(options)
    
    ## time stamp
    start_time_partial = HCGB_time.timestamp(start_time_total)
    
    ## debug message
    if (Debug):
        HCGB_aes.debug_message("retrieve_database")
        HCGB_main.print_all_pandaDF(retrieve_databases)
    
    ## init
    dataFrame_MLST = pd.DataFrame()
    
    print ("+ Generate an species typification for each sample retrieved using:")

   ###########################################################################
    if options.species2use or options.other_MLST_profile:
    ###########################################################################
     
        print ("+ No species typification will be generated, as user provided some information:")
        HCGB_aes.warning_message("Information provided will be applied for all samples")
           
        ## get MLST_profile: default or provided
        mlst_profile_dict = MLST_caller.get_MLST_profiles()
       
        #name, genus, species,
        dataFrame_MLST = pd.DataFrame(columns=("sample", "genus", "species", "mlst"))
        
        if options.species2use:
            
            ## get string
            print ("+ Species provided by user: " + options.species2use)
            species2use = options.species2use.split(" ")
            mlst2use = (species2use[0][0] + species2use[1]).lower()
            
            try:
                if mlst_profile_dict[mlst2use]:
                    print(colored("\t- Species name matches MLST available: OK", 'green'))
            except:
                print("\n")
                HCGB_aes.error_message("MLST provided is not available")
                ## exit if error
                HCGB_aes.raise_and_exit("Check your mlst configuration or options with --help_MLST")
            
            ## species provided by user as option
             
            # Group dataframe sample name
            sample_results = pd_samples_retrieved.groupby(["name_sample"])
            for name, grouped in sample_results:
                dataFrame_MLST.loc[len(dataFrame_MLST)] = (name[0], species2use[0], species2use[1], mlst2use)
            ##########################################################################
    
        ###########################################################################
        elif options.other_MLST_profile:
        ###########################################################################
            ## user provides the name of the MLST to use
            
            ## get string
            MLST_profile2use = options.other_MLST_profile
            
            print ("+ MLST provided by user: " + MLST_profile2use)
            ## species provided by user as option
            
            try:
                if mlst_profile_dict[MLST_profile2use]:
                    print(colored("\t- MLST provided is OK", 'green'))
            except:
                print("\n")
                HCGB_aes.error_message("MLST provided is not available")
                ## exit if error
                HCGB_aes.raise_and_exit("Check your mlst configuration or options with --help_MLST")
   
            # Group dataframe sample name
            sample_results = pd_samples_retrieved.groupby(["name_sample"])
            species2use_list = mlst_profile_dict[MLST_profile2use].split(" ")
            for name, grouped in sample_results:
                dataFrame_MLST.loc[len(dataFrame_MLST)] = (name[0], species2use_list[0], species2use_list[1], MLST_profile2use)
            ###########################################################################
    
    
        if Debug:
            HCGB_aes.debug_message("dataFrame_edirect for MLST call")
            HCGB_main.print_all_pandaDF(dataFrame_MLST)

    ###########################################################################
    ## Determine species using kraken2
    elif options.kraken2: 
    ###########################################################################

        ## Flag to identify from scratch using kraken2

        ######## Kraken2 identification
        dataFrame_kraken = Kraken_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, start_time_partial)
        
        ## functions.timestamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)
        
        ## debug message
        if (Debug):
            HCGB_aes.debug_message("retrieve results to summarize")
            HCGB_aes.debug_message("dataFrame_kraken")
            HCGB_main.print_all_pandaDF(dataFrame_kraken)
            
            
        ## For later MLST
        dataFrame_MLST = dataFrame_kraken
            
    ###########################################################################
    elif options.kma:
    ###########################################################################
        
        HCGB_aes.raise_and_exit("KMA identification process is not enabled due to unexpected errors. Further debugging required!")
        exit()
    
        ## ATTENTION: This option was producing errors, it might need debugging
        ## KMA is not working so far, no valid output results are generated and continuos memory crash errors are produced...
        ## just skip the process
    
        ## Flag to identify from scratch using KMA
        ## using other databases, user provided sequences or input genbank entries
    
        ######## KMA identification
        dataFrame_kma = KMA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, start_time_partial)
        
        ## functions.timestamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)
        
        ## debug message
        if (Debug):
            HCGB_aes.debug_message("retrieve results to summarize")
            HCGB_aes.debug_message("dataFrame_kma")
            HCGB_main.print_all_pandaDF(dataFrame_kma)
            
        ###########################################################################
        ## exit if viral search
        skip=False
        if (len(options.kma_dbs) == 1):
            for i in options.kma_dbs:
                if (i == 'viral'):
                    print ()
                    MLST_results = ''
                    options.slow = False
                    skip=True
                
                ## what if only plasmids?
        ###########################################################################
        
        ## create input edicrect df
        ## Call edirect to retrieve selected bacteria identified
        ## do edirect and MLST if bacteria
        if (not skip):        
            dataFrame_edirect = pd.DataFrame()
            
            ######## EDirect identification
            dataFrame_MLST = edirect_ident(dataFrame_kma, outdir_dict, Debug)
            
            ## functions.timestamp
            start_time_partial = HCGB_time.timestamp(start_time_partial)
        
            ## debug message
            if (Debug):
                HCGB_aes.debug_message("retrieve results from NCBI")
                HCGB_main.print_all_pandaDF(dataFrame_edirect)

    ###########################################################################


    ###########################################################################
    # MLST identification
    ###########################################################################
    MLST_results = MLST_ident(options, pd_samples_retrieved, dataFrame_MLST, outdir_dict, start_time_partial)

    ## functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)

    ## debug message
    if (Debug):
        HCGB_aes.debug_message("MLST_results")
        HCGB_main.print_all_pandaDF(MLST_results)


    ###########################################################################
    ## generate summary for sample: all databases
    ## MLST, plasmids, genome, etc
    ###########################################################################
    HCGB_aes.boxymcboxface("Results Summary")
    
    #####################################
    ## Summary identification results  ##
    #####################################
    ## parse results
    if options.project:
        final_dir = os.path.join(outdir, 'report', 'ident')
        HCGB_files.create_folder(final_dir) 
    else:
        final_dir = outdir

    ###
    excel_folder = HCGB_files.create_subfolder("samples", final_dir)
    print ('+ Print summary results in folder: ', final_dir)
    print ('+ Print sample results in folder: ', excel_folder)
    
    dataFrame_MLST = dataFrame_MLST.rename(columns={"new_name": "name_sample"})
    
    ## debug message
    if (Debug):
        print("Merge")
        HCGB_aes.debug_message("MLST_results")
        HCGB_main.print_all_pandaDF(MLST_results)
    
        HCGB_aes.debug_message("dataFrame_MLST")
        HCGB_main.print_all_pandaDF(dataFrame_MLST)

    # Merge dataframe results summary by sample name
    sample_results_summary = pd.merge(MLST_results, dataFrame_MLST, on='name_sample')
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: sample_results_summary merged **", 'yellow'))
        HCGB_main.print_all_pandaDF(sample_results_summary)
    
    ## init some variables
    bracken_results = pd.DataFrame()
    results_summary_KMA = pd.DataFrame()

    ###########################################################################
    ## Loop for each sample and save results in ident/samples folder generated
    ###########################################################################    
    for index, grouped in sample_results_summary.iterrows():
    ###########################################################################
    
        name = grouped['name']
        
        print("Printing information for sample: ")
        print(name)
    
        ## create a excel and txt for sample
        name_sample_excel = os.path.join(excel_folder,  name + '_ident.xlsx')
        name_sample_csv = os.path.join(outdir_dict[name], 'ident_summary.csv') ## check in detached mode

        ## A summary for all samples
        with pd.ExcelWriter(name_sample_excel, engine="xlsxwriter", engine_kwargs={"options": {"nan_inf_to_errors": True}}) as writer_sample:

            ###########################################################################
            if options.kma:
            ###########################################################################
                
                ## subset dataframe    & print result
                results_summary_toPrint_sample = grouped[['Sample','#Template',
                                                        'Query_Coverage','Template_Coverage',
                                                        'Depth', 'Database']] 
                results_summary_toPrint_sample.to_excel(writer_sample, sheet_name="KMA") ## write excel handle
                results_summary_toPrint_sample.to_csv(name_sample_csv) ## write csv for sample
                
            ###########################################################################
            elif options.kraken2:
            ###########################################################################
                
                bracken_df = pd.read_csv(grouped['sample'], sep="\t")
                bracken_df['name'] = grouped['name_sample']
                bracken_results = pd.concat([bracken_results, bracken_df])
    
                ## print kraken2/bracken results in xlsx format
                bracken_df.to_excel(writer_sample, sheet_name="Kraken2") ## write excel handle
                ## kraken already saves species name in .species file                
                
                
            ## read MLST
            if not MLST_results.empty:
                results_summary_MLST = grouped[['name_sample','species_id', 'scheme', 'alleles']] 
                results_summary_MLST.to_excel(writer_sample, sheet_name="MLST") ## write excel handle
    
    
    ###########################################################################
    ## A summary for all samples
    ###########################################################################
    name_excel = os.path.join(final_dir, 'identification_summary.xlsx')
    print ('+ Summary information in excel file: ', name_excel)
    with pd.ExcelWriter(name_excel, engine="xlsxwriter", engine_kwargs={"options": {"nan_inf_to_errors": True}}) as writer:
        ## write MLST
        if not MLST_results.empty:
            MLST_results.to_excel(writer, sheet_name='MLST')
    
        ## kma results
        if options.kma:            
            ## KMA dataframe: print result for sources
            results_summary_KMA = dataFrame_MLST[['Sample','#Template',
                                                'Query_Coverage','Template_Coverage',
                                                'Depth', 'Database']] 
            
            ## Sum plasmid and chromosome statistics ##
            ## sum coverage
            total_coverage = results_summary_KMA.groupby('Sample')['Query_Coverage'].sum().reset_index()
            
            ## debug message
            if (Debug):
                print ("*** Sum: Query_coverage ***")        
                print (total_coverage)
            
            ## TODO: Fix this chunk of code
            ## TODO: FIX SUMMARY REPORT
            results_summary_KMA = results_summary_KMA.set_index('Sample')
            results_summary_KMA = results_summary_KMA.sort_values(by=['Sample', 'Database', 'Query_Coverage'],ascending=[True, True,True])
            results_summary_KMA.to_excel(writer, sheet_name='KMA') ## write excel handle
        
            if (Debug):
                print (colored("**DEBUG: results_summary_KMA **", 'yellow'))
                HCGB_main.print_all_pandaDF(results_summary_KMA)

        ###
        elif options.kraken2:
            
            ## read ecah bracken file and create merge dataframe
            print()
            bracken_results.to_excel(writer, sheet_name='Kraken2')
            
            
        else:
            ## user provided species or mlst profile
            print()
            dataFrame_MLST.to_excel(writer, sheet_name='KMA') ## write excel handle

    ### timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)
    
    ######################################
    ## update database for later usage
    ######################################
    if options.slow:

        HCGB_aes.boxymcboxface("Update Sample Database")

        ## update db
        print ("+ Update database with samples identified")

        ## debug message
        if (Debug):
            print (colored("**DEBUG: dataFrame_edirect **", 'yellow'))
            pd.set_option('display.max_colwidth', None)
            pd.set_option('display.max_columns', None)
            print (dataFrame_edirect)

        ## dataFrame_edirect
        file_toprint = final_dir + '/edirect_info2download.csv'
        dataFrame_edirect.to_csv(file_toprint)

        ## update database with samples identified
        data2download = dataFrame_edirect.filter(['genus','species', 'strain', 'genome'])
        data2download = data2download.rename(columns={'genome': 'NCBI_assembly_ID', 'strain' : 'name'})
        NCBI_folder = os.path.abspath(options.database) + '/NCBI'
        database_generator.NCBI_DB(data2download, NCBI_folder, Debug)

    else:
        print ("+ No update of the database has been requested [Default]. Use option --slow instead")
        
    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## dump information and parameters
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"ident", "time":HCGB_time.timestamp(time.time()),
                "BacterialTyper version":pipeline_version }

    HCGB_info.dump_info_run(info_dir, 'ident', options, runInfo, options.debug)
    
    ## dump conda details
    HCGB_info.dump_info_conda(info_dir, "ident", options.debug)

    print ("+ Exiting identification module.")
    return()

####################################
def KMA_ident(options, pd_samples_retrieved, outdir_dict, retrieve_databases, time_partial):
    """Kmer identification using software KMA_.
    
    :param options: options passed to the :func:`BacterialTyper.modules.ident.run_ident` main function (threads, KMA_cutoff, etc). See details in...
    :param pd_samples_retrieved: pandas dataframe for samples to process.
    :param outdir_dict: dictionary containing information for each sample of the output folder for this process.
    :param retrieve_databases: 
    :param time_partial: timestamp of start time of the process.
    
    :type options: 
    :type pd_samples_retrieved: pandas.DataFrame()
    :type outdir_dict: Dictionary
    :type retrieve_databases: pandas.DataFrame()
    :type time_partial: 
    
    :return: Information of the identification. See example below.
    :rtype: pandas.DataFrame()
    
    See example of returned dataframe in file :file:`/devel/results/KMA_ident_example.csv` here:
    
    .. include:: ../../devel/results/KMA_ident_example.csv
        :literal:
    
    .. seealso:: This function depends on other ``BacterialTyper`` functions called:
    
        - :func:`BacterialTyper.config.set_config.get_exe`
    
        - :func:`BacterialTyper.scripts.functions.boxymcboxface`
        
        - :func:`BacterialTyper.modules.ident.send_kma_job`
        
        - :func:`BacterialTyper.modules.ident.get_outfile`
    
        - :func:`BacterialTyper.scripts.species_identification_KMA.check_db_indexed`
    
        - :func:`BacterialTyper.scripts.species_identification_KMA.parse_kma_results`
    
        
    .. include:: ../../links.inc    
    
    """    
    
    ### print header
    HCGB_aes.boxymcboxface("KMA Identification")

    ## set defaults
    kma_bin = set_config.get_exe("kma")    

    ## check status
    databases2use = []
    for    index, db2use in retrieve_databases.iterrows():
        ## index_name
        if (str(db2use['source']).startswith('KMA')):
            print ('+ Check database: ' + db2use['db'])
            fold_name = os.path.dirname(db2use['path'])
            
            index_status = KMA_caller.check_db_indexed(db2use['path'], fold_name )
            if (index_status == True):
                print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
                databases2use.append(db2use['path'])
            else:
                #databases2use.remove(db2use)
                print (colored("\t**Databases %s is not correctly indexed. Not using it...\n" % db2use['db'], 'red'))

    ## debug message
    if (Debug):
        print (colored("**DEBUG: databases2use\n" +  "\n".join(databases2use) + "\n**", 'yellow'))

    ## Start identification of samples
    print ("\n+ Send KMA identification jobs...")
    
    
    ### ATTENTION:
    ## KMA is not working so far, no valid output results are generated and continuos memory 
    ## crash errors are produced... 
    ## just skip the process
    results_summary = pd.DataFrame()
#    return results_summary

    ## optimize threads
    name_list = set(pd_samples_retrieved.index)
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["name_sample"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        for db2use in databases2use:
            try:
                print ("+ Removing database from memory...")
                return_code_rm = KMA_caller.remove_db(kma_bin, db2use)
            except:
                pass
            
            
            ## load database on memory
            print ("+ Loading database on memory for faster identification.")
            return_code_load = KMA_caller.load_db(kma_bin, db2use)
            
            print ("+ Sending jobs for species identification.")
            ## send for each sample
            commandsSent = { executor.submit(send_kma_job, 
                                            outdir_dict[name[0]], 
                                            sorted(cluster[ cluster['ext']=='fastq']["sample"].tolist()), 
                                            name[0], db2use, threads_job, Debug): name[0] for name, cluster in sample_frame }

            for cmd2 in concurrent.futures.as_completed(commandsSent):
                details = commandsSent[cmd2]
                try:
                    data = cmd2.result()
                except Exception as exc:
                    print ('***ERROR:')
                    print (cmd2)
                    print('%r generated an exception: %s' % (details, exc))

            ## remove database from memory
            print ("+ Removing database from memory...")
            return_code_rm = KMA_caller.remove_db(kma_bin, db2use)
            
            if (return_code_rm == 'FAIL'):
                cmd_rm_db = "kma shm -t_db %s -shmLvl 1 -destroy" %db2use
                print (colored("***ERROR: Removing database from memory failed. Please do it manually! Execute command: %s" %cmd_rm_db,'red'))
        
            ## functions.timestamp
            time_partial = HCGB_time.timestamp(time_partial)
        
    ## parse results        
    print ("+ KMA identification call finished for all samples...")
    print ("+ Parse results now")
    results_summary = pd.DataFrame()
    
    #return results_summary
    
    for db2use in databases2use:
        ### [TODO]: parse data according to database: bacteria, plasmids or user data or genbank data provided
        
        basename_db = os.path.basename(db2use)
        pd.set_option('display.max_colwidth', None)
        pd.set_option('display.max_columns', None)

        ###
        for name, cluster in sample_frame:
            
            ## get result
            ## outdir_KMA
            outdir_dict_kma = HCGB_files.create_subfolder("kma", outdir_dict[name])
            result = get_outfile(outdir_dict_kma, name, db2use)
            #print ('\t- File: ' + result + '.spa')
            
            ## get results using a cutoff value [Defaulta: 80]
            results = KMA_caller.parse_kma_results(result + '.spa', options.KMA_cutoff)
            results['Database'] = basename_db

            ### check if db2use is plasmids as it could be several.
            if (results.index.size > 1):
                if (basename_db == "plasmids.T" or basename_db == "viral.TG"):
                    ## let it be several entries
                    results['Sample'] = name
                    results_summary = results_summary.append(results, ignore_index=True)
                else:
                    print (colored("###########################################", 'yellow'))
                    print (colored("Sample %s contains multiple strains." %name, 'yellow'))
                    print (colored("###########################################", 'yellow'))
                    print (colored(results, 'yellow'))
                    print ('\n\n')
                    
                    ## add both strains if detected    
                    results['Sample'] = name
                    results_summary = results_summary.append(results, ignore_index=True)

                    ## TODO: Fix this chunk of code
                    ## TODO: add multi-isolate flag
        
            elif (results.index.size == 1): ## 1 clear reference
                results['Sample'] = name
                results_summary = results_summary.append(results, ignore_index=True)
        
            else:
                print (colored('\tNo clear strain from database %s has been assigned to sample %s' %(basename_db, name), 'yellow'))
                ## add empty line if no available
                results['Sample'] = name
                results_summary = results_summary.append(results, ignore_index=True)
    
    print ("+ Finish this step...")
    
    ## debug message
    if (Debug):
        results_summary.to_csv(quotechar='"')
    
    return (results_summary)

###################################
def send_kma_job(outdir_file, list_files, name, database, threads, Debug):
    """Executes KMA identification jobs
    
    This function automates the process of checking if any previous run succeeded or
    runs the appropiate identification process for the sample and database provided.
    
    :param outdir_file:
    :param list_files:
    :param name:
    :param database:
    :param threads:
    :param dataFrame_sample:
    
    :type outdir_file:
    :type list_files:
    :type name:
    :type database:
    :type threads:
    :type dataFrame_sample:
    
    .. seealso:: This function depends on other ``BacterialTyper`` functions called:
    
        - :func:`BacterialTyper.config.set_config.get_exe`
    
        - :func:`BacterialTyper.scripts.species_identification_KMA.kma_ident_call`
    
        - :func:`BacterialTyper.module.ident.get_outfile`
        
        - :func:`BacterialTyper.scripts.functions.read_time_stamp`
        
        
    """
    
    if (Debug):
        print (colored("**DEBUG: ident.send_kma_job call**", 'yellow'))
        print ("outdir_file")
        print (outdir_file)
        print ("list_files")
        print (list_files)
        print ("name: " + name)
        print ("database: " + database)
        
    ## outdir_KMA
    outdir_dict_kma = HCGB_files.create_subfolder("kma", outdir_file)

    ## set defaults
    kma_bin = set_config.get_exe("kma")

    ## get outfile
    outfile = get_outfile(outdir_dict_kma, name, database)

    ## check if previously run and succeeded
    basename_tag = os.path.basename(outfile)
    filename_stamp = outdir_dict_kma + '/.success_' + basename_tag
    
    if (Debug):
        print ("Outdir: ", outdir_dict_kma)
        print ("outfile: ", outfile)
        print ("Filename_stamp: ", filename_stamp)
    
    if os.path.isfile(filename_stamp):
        stamp =    HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
    else:
        ## debug message
        if (Debug):
            print (colored("**DEBUG: species_identification_KMA.kma_ident_module call**", 'yellow'))
            print ("outfile = get_outfile(outdir_dict_kma, name, db2use)")
            print ("outfile: ", outfile)
            print ("species_identification_KMA.kma_ident_module(outfile, list_files, name, database, threads) ")
            print ("species_identification_KMA.kma_ident_module" + "\t" + outfile + "\t" + str(list_files) + "\t" + name + "\t" + database + "\t" + str(threads) + "\n") 
    
        ## Sparse or not
        #if any(name in basename_tag for name in ['userData_KMA', 'genbank_KMA']):
#        if (basename_tag == 'userData_KMA'):
#            option = ''
#        else:
#            option = '-Sparse '
    
        ## Add option to retrieve databse from memory
        option = ""
        option = option + '-shm 1'
    
        # Call KMA
        KMA_caller.kma_ident_call(outfile, list_files, name, database, kma_bin, option, threads)
        stamp =    HCGB_time.print_time_stamp(filename_stamp)

####################################
def get_outfile(output_dir, name, index_name):
    """
    Generates the name for the output file created
    
    :param output_dir: Absolute path to results folder 
    :param name: Name of the sample
    :param index_name: Name of the database
    
    :type output_dir: string
    :type name: string
    :type index_name: string
    
    :retruns: Output file absolute path
    """
    basename_tag = os.path.basename(index_name)
    if Project:
        output_path = output_dir
    else:
        output_path = HCGB_files.create_subfolder(name, output_dir)
        
    out_file = output_path + '/' + name + '_' + basename_tag    
    return(out_file)
   
####################################
def edirect_ident(dataFrame, outdir_dict, Debug):
    """Connect to NCBI for information retrieval
    
    This functions uses the software edirect_ to connect to NCBI and retrieve some information regarding samples, assemblies, publications, etc.
    
    :param dataFrame: pandas dataframe for samples to process. Result from :func:`BacterialTyper.modules.ident.KMA_ident`.
    :param outdir_dict: dictionary containing information for each sample of the output folder for this process.
    
    :type dataFrame: pandas.DataFrame()
    :type outdir_dict: Dictionary
    
    :return: Information of the identification 
    :rtype: pandas.DataFrame()
    
    See example of returned dataframe in file :file:`/devel/results/edirect_download_results.csv` here:
    
    .. include:: ../../devel/results/edirect_download_results.csv
        :literal:
    
    .. seealso:: This function depends on other ``BacterialTyper`` functions called:
    
        - :func:`BacterialTyper.scripts.functions.get_info_file`
        
        - :func:`BacterialTyper.scripts.functions.read_time_stamp`
    
        - :func:`BacterialTyper.scripts.functions.print_time_stamp`

        - :func:`BacterialTyper.scripts.functions.optimize_threads`
    
        - :func:`BacterialTyper.scripts.functions.create_subfolder`
    
        - :func:`BacterialTyper.scripts.functions.boxymcboxface`
        
        - :func:`BacterialTyper.scripts.functions.is_non_zero_file`
    
        - :func:`BacterialTyper.scripts.edirect_caller.generate_docsum_call`
        
        - :func:`BacterialTyper.scripts.edirect_caller.generate_xtract_call`
        
    .. include:: ../../links.inc    
    """
    
    
    ### ATTENTION:
    ## KMA is not working so far, no valid output results are generated and continuos memory crash errors are produced...
    ## just skip the process
    return pd.DataFrame()

    
    ################################################
    ## TODO: Fix this chunk of code
    ## TODO: What to do if multi-isolate sample?
    ################################################
    
    ## edirect    
    HCGB_aes.boxymcboxface("EDirect information")
    print ("+ Connect to NCBI to get information from samples identified...")

    ## create dataframe to return results
    edirect_frame = pd.DataFrame(columns=("sample", "genus", "species", "strain", "BioSample", "genome", "Plasmids"))

    ## debugging messages
    if Debug:
        print ("*******************************************************")
        print ("Dataframe sample_results: ")
        
    # Group dataframe sample name
    sample_results = dataFrame.groupby(["Sample"])
    
    for name, grouped in sample_results:
        ## debugging messages
        if Debug:
            print ("Name: ", name)
            print (grouped)
        
        ## use edirect to get Species_name and entry for later identification
        edirect_folder = HCGB_files.create_subfolder('edirect', outdir_dict[name])
        
        ## chromosome match
        if (len(grouped.loc[grouped['Database'] == 'bacteria.ATG']['#Template']) == 0):
            if Debug:
                print ("Name: ", name)
                print ("No chromosome match identified by kmer")
            
            genus = ''
            species = ''
            BioSample_name = ''
            AssemblyAcc = ''            
        
        else:
            nucc_entry = grouped.loc[grouped['Database'] == 'bacteria.ATG']['#Template'].values[0].split()
            ## e.g. NZ_CP029680.1 Staphylococcus aureus strain AR_0215 chromosome, complete genome

            ##
            out_docsum_file = edirect_folder + '/nuccore_docsum.txt'
            tmp_species_outfile = edirect_folder + '/info.csv'
            filename_stamp = edirect_folder + '/.success_species'
                    
            if os.path.isfile(filename_stamp):
                stamp =    HCGB_time.read_time_stamp(filename_stamp)
                print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
                status=True    
            else: 
                edirect_caller.generate_docsum_call('nuccore', nucc_entry[0], out_docsum_file)
                status = edirect_caller.generate_xtract_call(out_docsum_file, 'DocumentSummary', 'Organism,BioSample,AssemblyAcc,Strain', tmp_species_outfile)
                
            ########################################    
            ## get information from edirect call
            ########################################    
            if not status:
                print ("NO INFORMATION")
                continue

            taxa_name_tmp = HCGB_main.get_info_file(tmp_species_outfile)
            Organism = taxa_name_tmp[0].split(',')[0].split()
            genus = Organism[0]                             ## genus
            species = Organism[1]                             ## species
            BioSample_name = taxa_name_tmp[0].split(',')[1]    ## BioSample
            AssemblyAcc = taxa_name_tmp[0].split(',')[2]     ## AssemblyAcc
            
            ## sometimes strain is missing
            if len(taxa_name_tmp[0].split(',')) > 3:
                strain = taxa_name_tmp[0].split(',')[3]     ## strain
            else:
                strain = 'NaN'
            
            ## get GenBank accession ID
            out_docsum_file_assembly = edirect_folder + '/assembly_docsum.txt'
            AssemblyAcc_outfile = edirect_folder + '/AssemblyAcc.csv'
            
            edirect_caller.generate_docsum_call('assembly', AssemblyAcc, out_docsum_file_assembly)
            edirect_caller.generate_xtract_call(out_docsum_file_assembly, 'DocumentSummary', 'Genbank', AssemblyAcc_outfile) 
            
            ## some error occurred
            if not HCGB_main.is_non_zero_file(out_docsum_file_assembly):
                continue
            
            ## Is it better to download Refseq or Genbank?
            ## https://www.quora.com/What-is-the-difference-between-Refseq-and-Genbank        
            
            GenbankAcc = HCGB_main.get_info_file(AssemblyAcc_outfile)
            if Debug:
                print("Sample: ", name)
                print("Genbank Acc: ", GenbankAcc[0])
        
        ## plasmid match
        group_plasmid = grouped.loc[grouped['Database'] == 'plasmids.T' ]
        plasmid_entries = group_plasmid['#Template'].tolist()
            ## e.g. NZ_CP029083.1 Staphylococcus aureus strain AR464 plasmid unnamed1, complete sequence
        plasmid_entries_str = ",".join([i.split()[0] for i in plasmid_entries])

        ## save edirect_frame
        #("sample", "taxa", strain, genome "BioSample", "Plasmids"))
        edirect_frame.loc[len(edirect_frame)] = (name, genus, species, strain, BioSample_name, GenbankAcc[0], plasmid_entries_str)

        stamp =    HCGB_time.print_time_stamp(filename_stamp)

    ## debugging messages
    if Debug:
        print ("*******************************************************")
    
    return (edirect_frame)

####################################
def MLST_ident(options, dataFrame, dataFrame_species, outdir_dict, time_partial):
    """
    
    conda_env/db/pubmlst
    
    https://rest.pubmlst.org/db/...

    Parameters
    ----------
    options : TYPE
        DESCRIPTION.
    dataFrame : TYPE
        DESCRIPTION.
    outdir_dict : TYPE
        DESCRIPTION.
    dataFrame_edirect : TYPE
        DESCRIPTION.
    retrieve_databases : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    ## MLST call    
    HCGB_aes.boxymcboxface("MLST typing")
    print ("+ Create classical MLST typification of each sample according to species retrieved by kmer...")
    
    ## get assembly files
    subset_Df = dataFrame[ dataFrame['tag'] == 'assembly']
    
    ## get MLST_profile: default or provided
    mlst_profile_dict = MLST_caller.get_MLST_profiles()
    
    ## debug message
    if (Debug):
        HCGB_aes.debug_message("dataFrame_species identified")
        print (dataFrame_species)
    
        HCGB_aes.debug_message("subset_Df")
        print (subset_Df)    
    
        HCGB_aes.debug_message("mlst_profile_dict")
        pprint.pprint(mlst_profile_dict)    


    ## Start identification of samples
    print ("\n+ Send MLST identification jobs...")
    
    ## optimize threads
    name_list = set(subset_Df.index)
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)
     
    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))
    
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
         ## send for each sample
         commandsSent = { executor.submit(MLST_caller.MLST_call,
                                          outfolder = outdir_dict[row['name_sample']],  ## outfolder
                                          assembly_file = row['sample'],             ## assembly file
                                          mlst_profile = dataFrame_species.loc[ dataFrame_species['name'] == row['name_sample'], 'mlst'].item(), 
                                          sample_name = row['name_sample'],
                                          minid=options.minid, 
                                          mincov=options.mincov, 
                                          minscore=options.minscore, 
                                          debug=Debug): index for index, row in subset_Df.iterrows() }
    
         for cmd2 in concurrent.futures.as_completed(commandsSent):
             details = commandsSent[cmd2]
             try:
                 data = cmd2.result()
             except Exception as exc:
                 print ('***ERROR:')
                 print (cmd2)
                 print(details)
                 print('%r generated an exception: %s' % (details, exc))
    
         ## functions.timestamp
         time_partial = HCGB_time.timestamp(time_partial)


    ## parse results and return
    print("+ Parse MLST results now...")
    MLST_results = pd.DataFrame()
    for sample_id in outdir_dict:
        print("\t+ Parsing: " + sample_id)
        
        json_file = os.path.join(outdir_dict[sample_id], 'MLST/MLST_res.json')
        sample_dict = HCGB_main.read_json_file(json_file, debug=Debug)
        if Debug:
            HCGB_aes.debug_message("Sample dictionary from json file:")
            print("JSON file: " + json_file)
            print("Information:")            
            print(sample_dict)
        
        ## concat
        MLST_results = pd.concat([MLST_results, pd.DataFrame.from_dict(sample_dict)])
        
    ## rename and return
    MLST_results.rename(columns = {'id':'name_sample'}, inplace=True)
    return(MLST_results)

####################################
def get_external_kma(kma_external_files, Debug):
    print ('\t- Get additional kma databases:')
    ## external sequences provided are indexed and generated in the same folder provided 
    
    option_db = ""
    if (kma_external_files):
        kma_external_files = set(kma_external_files)        
        kma_external_files = [os.path.abspath(f) for f in kma_external_files]    
        
        ## check if indexed and/or index if necessary
        external_kma_dbs_list = []
        
        ## set defaults
        kma_bin = set_config.get_exe("kma")
        for f in kma_external_files:
            file_name = os.path.basename(f)
            fold_name = os.path.dirname(f)
            print (colored('\t\t+ %s' %file_name, 'green'))
            print ()

            ## generate db
            databaseKMA = KMA_caller.generate_db([f], file_name, fold_name, 'new', 'single', Debug, kma_bin)
            if not databaseKMA:
                print (colored("***ERROR: Database provided is not indexed.\n" %databaseKMA,'orange'))
            else:
                external_kma_dbs_list.append(databaseKMA)
            
        external_kma_dbs_string = ','.join(external_kma_dbs_list)
        option_db = "kma_external:" + external_kma_dbs_string

    else:
        ## rise error & exit
        print (colored("***ERROR: No database provided via --kma_external_file option.\n",'red'))
        exit()

    return(option_db)                

####################################
def get_options_db(options):
    """Select databases to use according to the input options.
    
    :param options:
    
    :returns: Dataframe with database information among all databases available.
    """
    
    print ("\n\n+ Select databases to use for identification:")
    
    ### database folder to use
    database2use = os.path.abspath(options.database)
    
    ## debug message
    if (Debug):
        HCGB_aes.debug_message("Database to use: " +  database2use)
    
    ## according to user input: select databases to use
    option_db = ""
    
    ############################################################
    ## Default db KMA
    ############################################################
    kma_dbs = []
    if not options.only_kma_db: ## exclusive
        #kma_dbs = ["bacteria", "plasmids"]
        kma_dbs = ["bacteria"]
            
    if (options.kma_dbs):
        options.kma_dbs = options.kma_dbs + kma_dbs
        options.kma_dbs = set(options.kma_dbs)        
    else:
        options.kma_dbs = kma_dbs

    ## rise error & exit if no dbs provided
    if not (options.kma_dbs):
        HCGB_aes.raise_and_exit("No database provided via --kma_db option.\n")

    ############################################################
    ### Options:
    
    ############
    ## 1) only user data: previously identified and added
    ############
    if (options.only_user_data):
        option_db = "user_data"
        
    ############
    ## 2) only genbank data: previously download from NCBI reference genomes
    ############
    elif (options.only_genbank_data):
        option_db = "genbank"
    
    ############
    ## 3) only external kma
    ############
    elif (options.only_external_kma):
        option_db = database_generator.get_external_kma(options.kma_external_files, Debug)
        ## rise attention
        if (options.kma_dbs):
            HCGB_aes.warning_message("Defatult databases and databases provided via --kma_dbs \
                                     option would not be used as --only_external_kma option provided.\n")

    #####################
    ## all databases KMA
    #####################
    else:        
        ####################
        ## default KMA dbs
        ####################
        print ('\t- Selecting kma databases:')
        kma_dbs_string = ','.join(options.kma_dbs)
        option_db = "kma:" + kma_dbs_string
    
        for i in options.kma_dbs:
            print (colored('\t\t+ %s' %i, 'green'))
        
        #################
        ## External file
        #################
        if (options.kma_external_files):
            option_db_tmp = database_generator.get_external_kma(options.kma_external_files, Debug)
            option_db = option_db + '#' + option_db_tmp
            
        #############################
        ## Previously identified data
        #############################
        if any([options.user_data, options.all_data]):
            option_db = option_db + '#kma_user_data:user_data'
    
        #############################
        ## Genbank reference data
        #############################
        if any([options.genbank_data, options.all_data]):
            option_db = option_db + '#kma_NCBI:genbank'


    ##############
    ### get dbs
    ###############
    print ("\n+ Parsing information to retrieve databases")
    print ("+ Reading from database: " + database2use)
    HCGB_aes.print_sepLine("-",50, False)

    ## get KMA databases    
    pd_KMA = database_generator.getdbs("KMA", database2use, option_db, Debug)    

    #####################
    ## databases Kraken
    #####################
    if options.kraken_dbs:
        k2_db = os.path.join(options.database, "Kraken", options.kraken_dbs)
        if os.path.isdir(k2_db):
            option_db = option_db + '#kraken:' + options.kraken_dbs
        else:
            try:
                kraken2_caller.get_dbs(db_name=options.kraken_dbs, fold=os.path.join(options.database, "Kraken"))
            except:
                HCGB_aes.raise_and_exit("No database available provided for Kraken2")
        
    elif options.other_kraken2_db:
        option_db = option_db + '#user_kraken:' + os.path.abspath(options.other_kraken2_db)

    ###############
    ## debug message
    if (Debug):
        HCGB_aes.debug_message("option_db: " +  option_db)
        
    ## Get kraken databases
    pd_Kraken = database_generator.getdbs("Kraken2", database2use, option_db, Debug)    
    
    if (Debug):
        HCGB_aes.debug_message("pd_KMA")
        HCGB_main.print_all_pandaDF(pd_KMA)
        HCGB_aes.debug_message("pd_Kraken")
        HCGB_main.print_all_pandaDF(pd_Kraken)

    HCGB_aes.print_sepLine("-",50, False)

    return (pd.concat([pd_KMA, pd_Kraken]))

####################################
def Kraken_ident(options, dataFrame_samples, outdir_dict, retrieve_databases, time_partial):
    """
    

    Parameters
    ----------
    options : TYPE
        DESCRIPTION.
    dataFrame_samples : TYPE
        DESCRIPTION.
    outdir_dict : TYPE
        DESCRIPTION.
    retrieve_databases : TYPE
        DESCRIPTION.
    time_partial : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    ### print header
    HCGB_aes.boxymcboxface("Kraken2 Identification")

    ## check status
    databases2use = []
    for index, db2use in retrieve_databases.iterrows():
        ## index_name
        if (str(db2use['source']).startswith('kraken')):
            print ('+ Check database: ' + db2use['db'])
            print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
            databases2use.append(db2use['path'])
            
            #fold_name = os.path.dirname(db2use['path'])
            # index_status = KMA_caller.check_db_indexed(db2use['path'], fold_name )
            # if (index_status == True):
            #     print (colored("\t+ Databases %s seems to be fine...\n\n" % db2use['db'], 'green'))
            #     
            # else:
            #     #databases2use.remove(db2use)
            #     print (colored("\t**Databases %s is not correctly indexed. Not using it...\n" % db2use['db'], 'red'))

    ## debug message
    if (Debug):
        print (colored("**DEBUG: databases2use\n" +  "\n".join(databases2use) + "\n**", 'yellow'))

    ## Start identification of samples
    print ("\n+ Send Kraken2 identification jobs...")
    
    ## optimize threads
    name_list = set(dataFrame_samples.index)
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    # Group dataframe by sample name
    sample_frame = dataFrame_samples.groupby(["name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
         ## send for each sample
         ## sample_name, read_files, threads_num, db_fold, outfolder, 
         ## level_abundance="S", thres_count = 50, hit_groups=3, others="", debug=False
         
         commandsSent = { executor.submit(kraken2_caller.kraken_run_all,
                         name[0],                                                       ## sample_name
                         sorted(cluster[ cluster['ext']=='fastq']["sample"].tolist()),  ## read_files
                         threads_job,                                                   ## threads                    
                         databases2use[0],                                              ## database
                         outdir_dict[name[0]],                                          ## outfolder
                         options.kraken2_level,                                         ## level_abundance
                         options.kraken2_read_count,                                    ## read threshold
                         options.kraken2_hit_groups,                                    ## hit groups kraken2
                         options.others_kraken2,                                        ## other options
                         Debug): name[0] for name, cluster in sample_frame }

         for cmd2 in concurrent.futures.as_completed(commandsSent):
             details = commandsSent[cmd2]
             try:
                 data = cmd2.result()
             except Exception as exc:
                 print ('***ERROR:')
                 print (cmd2)
                 print('%r generated an exception: %s' % (details, exc))

         ## functions.timestamp
         time_partial = HCGB_time.timestamp(time_partial)
        
    ## parse results        
    print ("+ Kraken2 identification call finished for all samples...")
    print ("+ Parse results now")
    
    ## parse results and generate summary file
    results_summary = pd.DataFrame()
    pd_samples_bracken = sampleParser.files.get_files(options, options.input, "ident", ["bracken"], options.debug)

    ## debug message
    if (Debug):
        HCGB_aes.debug_message("pd_samples_bracken")
        HCGB_main.print_all_pandaDF(pd_samples_bracken)


    ## 
    mlst_profile_dict = MLST_caller.get_MLST_profiles()
    pd_samples_bracken['species_id'] = ""
    pd_samples_bracken['mlst'] = ""

    print("+ Parsing kraken2/bracken results")
    for i in range( len(pd_samples_bracken) ):
        print("Sample: " + pd_samples_bracken.loc[i, "name"])
        try:
            species2use = kraken2_caller.parse_results(bracken_res= pd_samples_bracken.loc[i, "sample"], 
                                     folder = os.path.dirname(pd_samples_bracken.loc[i, "sample"]), cut_off=0.90)
        except:
            HCGB_aes.error_message("Parsing error for Kraken output")
            print(pd_samples_bracken.loc[i])
            HCGB_aes.raise_and_exit("")

        print("\tSpecies: " + species2use)

        pd_samples_bracken.loc[i, "species_id"] = species2use
        species2use_list = species2use.split(" ")
        mlst2use = (species2use_list[0][0] + species2use_list[1]).lower()
        print("\tMLST profile: " + mlst2use)
        try:
            if mlst_profile_dict[mlst2use]:
                print(colored("\t- Species name matches MLST available: OK", 'green'))
                pd_samples_bracken.loc[i, "mlst"] = mlst2use

        except:
            print("\n")
            HCGB_aes.error_message("MLST provided is not available")
            HCGB_aes.warning_message("This sample will not be processed in the MLST analysis")

    ## debug message
    if (Debug):
        HCGB_aes.debug_message("pd_samples_bracken")
        HCGB_main.print_all_pandaDF(pd_samples_bracken)
    
    return pd_samples_bracken
    
