###############
def BUSCO_download(name, ftp, folder):
    """Downloads BUSCO datasets provided
    
    This code checks if dataset is already available in folder provided. If not available proceeds to download it.
    
    It creates a subfolder (using :func:`BacterialTyper.scripts.functions.create_folder`) and downloads the dataset from ftp site provided (using :func:`BacterialTyper.scripts.functions.wget_download`).

    The file donwloaded would gunzipped so it is decompressed (using  :func:`BacterialTyper.scripts.functions.extract`). 
    
    A timestamp is printed to reflect download data (using :func:`BacterialTyper.scripts.functions.print_time_stamp`).
    
    The data download is checked for integrity using the function :func:`BacterialTyper.scripts.BUSCO_caller.BUSCO_check_dataset`.

    If the process failed, the whole process is retried (CAUTION: it might generate infinite loop).    
    
    :param name: Dataset name provided.
    :param ftp: FTP site for dataset.
    :param folder: Absolute path to folder to store results (e.g. database/BUSCO).
    :type name: string
    :type ftp: string
    :type folder: string
    
    :returns: Folder absolute path containing downloaded dataset.
    :rtype: string
    :warnings: Returns **FAIL** if check process failed.
    
    
    .. seealso:: This function depends on other BacterialTyper functions called:
    
        - :func:`BacterialTyper.scripts.functions.create_folder`
        
        - :func:`BacterialTyper.scripts.functions.wget_download`
        
        - :func:`BacterialTyper.scripts.functions.extract`
        
        - :func:`BacterialTyper.scripts.functions.print_time_stamp`
        
        - :func:`BacterialTyper.scripts.functions.print_sepLine`

        - :func:`BacterialTyper.scripts.BUSCO_caller.BUSCO_check_dataset`

    """
    
    print (colored("\n+ BUSCO dataset: " + name + " - v9 OrthoDB", 'yellow')) ## modify OrthoDB version if changed
    subfolder = folder + '/' + name
    file_name = os.path.basename(ftp)
    path_file = subfolder + '/' + file_name
    folderName = subfolder + '/' + file_name.split('.tar.gz')[0]

    if os.path.exists(subfolder):
        print ('Already available in the path provided...') 
    else:
        ## does not exists
        functions.create_folder(subfolder)

        functions.wget_download(ftp, subfolder)

        ## extract
        print ("+ Extract file...")
        functions.extract(path_file, subfolder)

        ## timestamp
        filename_stamp = subfolder + '/.success'
        functions.print_time_stamp(filename_stamp)

    ## check if it is allright
    code = BUSCO_check_dataset(folderName)
    functions.print_sepLine("-", 50, False)
    
    if (code == 'FAIL'):
        print (colored('*** Dataset failed. Try to download it again...','red'))
        shutil.rmtree(folder)
        BUSCO_download(name, ftp, folder)
        
        ### CAUTION: check for infinite loop
        
    return (folderName)

