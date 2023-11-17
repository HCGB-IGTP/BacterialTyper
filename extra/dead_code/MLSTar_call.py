#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:30:17 2023

@author: jsanchez
"""

####################################
def MLSTr_ident(options, dataFrame, outdir_dict, dataFrame_edirect, retrieve_databases):
    """Generate MLST profile identification
    
    This functions uses the `MLSTar software`_ to retrieve Multi locus sequence typing (MLST) profiles from PubMLST_ for the given species previously identified by KMA. It generates MLST profiling for each sample. 
    
    :param options: options passed to the :func:`BacterialTyper.modules.ident.run_ident` main function (threads, KMA_cutoff, etc). See details in...
    :param dataFrame: pandas dataframe for samples to process. Result from :func:`BacterialTyper.modules.ident.KMA_ident`.
    :param outdir_dict: dictionary containing information for each sample of the output folder for this process.
    :param dataFrame_edirect: pandas dataframe resulted from :func:`BacterialTyper.modules.ident.edirect_ident`.
    :param retrieve_databases: 
    
    :type options: 
    :type dataFrame: pandas.DataFrame()
    :type outdir_dict: Dictionary
    :type dataFrame_edirect: pandas.DataFrame()
    :type retrieve_databases: pandas.DataFrame()
    
    :return: Information of the MLST identification. Dictionary keys are samples and values are the absolute path to file generate by :func:`BacterialTyper.scripts.MLSTar.run_doMLST` containing MLST information.
    :rtype: Dictionary

    
    See example of returned dataframe in file :file:`/devel/results/doMLST_result_example.csv` here:
    
    .. include:: ../../devel/results/doMLST_result_example.csv
        :literal:
    
    .. seealso:: Additional information to PubMLST available datasets.
    
        - :doc:`PubMLST datasets<../../../data/PubMLST_datasets>`
    
    
    .. seealso:: This function depends on other ``BacterialTyper`` functions called:
    
        - :func:`BacterialTyper.scripts.functions.read_time_stamp`
    
        - :func:`BacterialTyper.scripts.functions.create_subfolder`
        
        - :func:`BacterialTyper.scripts.functions.boxymcboxface`
        
        - :func:`BacterialTyper.scripts.MLSTar.run_MLSTar`
        
        - :func:`HCGB.sampleParser.files.get_files`
        
        - :func:`BacterialTyper.scripts.MLSTar.get_MLSTar_species`
        
    .. include:: ../../links.inc    
    """
    ## set config
    rscript = set_config.get_exe("Rscript")
    
    ## TODO: Samples might not be assembled...to take into account and return 0
    
    ## TODO: Fix and install MLSTar during installation
    
    print(colored("Attention: Fix:", 'red'))
    print(MLSTar.get_MLSTar_package_installed())
    
    ########################################################################################

    ## TODO: Fix this chunk of code
    ## TODO: What to do if multi-isolate sample?
    ## TODO: Control if a different profile is provided via --MLST_profile
    ## TODO: Check time passed and download again if >?? days passed]
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: dataFrame_edirect identified**", 'yellow'))
        print (dataFrame_edirect)

    ## MLST call    
    HCGB_aes.boxymcboxface("MLST typing")
    print ("+ Create classical MLST typification of each sample according to species retrieved by kmer...")

    ## get assembly files
    input_dir = os.path.abspath(options.input)
    assembly_samples_retrieved = sampleParser.files.get_files(options, input_dir, "assembly", ["fna"], options.debug)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: assembly_samples_retrieved**", 'yellow'))
        print (assembly_samples_retrieved)    
    
    # init
    MLST_results = {}
        
    ## get MLST_profile: default or provided
    mlst_profile_list = retrieve_databases.loc[ retrieve_databases['db'] == 'PubMLST']['path'].tolist()

    if (Debug):
        print ("** Debug **")
        print ("mlst_profile_list")
        print (mlst_profile_list)

        print ("dataFrame_edirect")
        print (dataFrame_edirect)


    ## Generate MLST call according to species identified for each sample
    for index, row in dataFrame_edirect.iterrows():
        MLSTar_taxa_name = MLSTar.get_MLSTar_species(row['genus'], row['species'] )
        
        if (MLSTar_taxa_name == 'NaN'):
            print (colored("\t- Not available PubMLST profile for sample [%s] identified as %s %s" %(row['sample'], row['genus'], row['species']), 'yellow'))
        
        else:
            for mlst_profile in mlst_profile_list:

                ## species folder
                #species_mlst_folder = functions.create_subfolder(MLSTar_taxa_name, pubmlst_folder)
                species_mlst = mlst_profile.split(',')[0]
                species_mlst_folder = mlst_profile.split(',')[1]
            
                ## output file
                output_file = species_mlst_folder + '/PubMLST_available_scheme.csv'
                filename_stamp = species_mlst_folder + '/.success_scheme'
            
                ## 
                if MLSTar_taxa_name == species_mlst:        
                    if os.path.isfile(filename_stamp):
                        stamp =    HCGB_time.read_time_stamp(filename_stamp)
                        print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
                    else:
                        ### get scheme available
                        MLSTar.getPUBMLST(MLSTar_taxa_name, rscript, output_file)
                        stamp =    HCGB_time.print_time_stamp(filename_stamp)

                    ## parse and get scheme for classical MLST
                    schemes_MLST = pd.read_csv(output_file, sep=',', header=0)

                    ##
                    for item, cluster in schemes_MLST.iterrows():
                        if cluster['len'] < 10:
                            scheme2use = int(cluster['scheme'])
                            continue            
                    ### 
                    sample = row['sample']
                    MLSTar_folder = HCGB_files.create_subfolder('MLST', outdir_dict[sample])
                    genome_file = assembly_samples_retrieved.loc[assembly_samples_retrieved['name'] == sample]['sample'].values[0]
    
                    ## call MLST
                    (results, profile_folder) = MLSTar.run_MLSTar(species_mlst_folder, rscript, MLSTar_taxa_name, scheme2use, sample, MLSTar_folder, genome_file, options.threads)
                    MLST_results[sample] = results

    ##    
    print ("+ Finish this step...")
    return (MLST_results)

