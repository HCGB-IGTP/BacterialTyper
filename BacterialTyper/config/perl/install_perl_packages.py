##################    
def perl_package_install_targz(name, version2install, http_tar_gz, install_dir, Debug):
    """
    Retrieves information for perl package installation
    
    This function installs packages archived in tar.gz files. It downloads and extracts
    ``tar.gz`` files and install them locally. This option might not be the best option
    as long as it would not install all sub-dependencies for each package. We would use
    preferentially :func:`BacterialTyper.config.install_dependencies.perl_package_install_cpan`.  

    :param name: Perl package name
    :param version2install: Version to install
    :param http_tar_gz: FTP/https site of the tar gz perl package (cpan)
    :param install_dir: Installation directory
    :param Debug: True/False for debugging messages

    :type name: string
    :type version2install: string 
    :type http_tar_gz: string 
    :type install_dir: string 
    :type Debug: boolean

    .. seealso: This function depends on other ``BacterialTyper`` functions such as:

        - :func:`BacterialTyper.scripts.functions.wget_download`

        - :func:`BacterialTyper.scripts.functions.extract`

        - :func:`BacterialTyper.config.install_dependencies.install_package`

    """    

    print (colored("Install missing perl package: " + name, 'yellow'))

    ## create folders
    perlPackages = functions.create_subfolder("perl_packages", install_dir)
    path2download_name = functions.create_subfolder(name, perlPackages)

    ## debugging messages
    if (Debug):
        print ("perlPackages: " + perlPackages)
        print ("path2download_name: " + path2download_name)

        ## download
    functions.wget_download(http_tar_gz, path2download_name)

    ## get targz file
    tar_gz = functions.retrieve_matching_files(path2download_name, 'tar.gz', Debug)
    functions.extract(tar_gz[0], path2download_name)

    ## debugging messages
    if (Debug):
        print ("** DEBUG: **")
        print ("http: " + http_tar_gz)
        print ('tar_gz: ')
        print (tar_gz)
        print ("*******************")

    ## install
    path2download_extract_folder = os.path.join(path2download_name, name)
    blib_lib = install_package_targz(path2download_extract_folder, install_dir, Debug, name)

    ## include in $PERL5LIB system variable
    print (colored("** ATTENTION:", 'yellow'))
    print ("Include the following path within your PERL5LIB variable:")
    print (blib_lib)
    
    ## debugging messages
    if (Debug):
        print ("** DEBUG: **")
        print ("blib_lib: " + blib_lib)

    return(version2install)

##################

#######################
def install_package_targz(package_path, install_path, Debug, name):
    """
    Install perl package downloaded
    
    :param package_path: Path where the perl package is extracted
    :param install_path: Path to install package

    :type package_path: string
    :type install_path: string

    :returns: Path where the module is installed to include in $PERL5LIB
    """
    ## change dir to package path
    os.chdir(package_path)

    print ("## Installing module: " + name + " ##")

    ## perl Makefile.PL
    makefile_perl = functions.retrieve_matching_files(package_path, "Makefile.PL", Debug)
    perl_exe = set_config.get_exe("perl", Debug)

    ## two installation possibilities: Makefile.PL or Build.PL
    if (makefile_perl):

        print ("+ Create make file")
        perl_MakeFile_cmd = perl_exe + ' ' + makefile_perl[0]

        ## debug messages
        if (Debug):
            print ("** Debug: chdir " + package_path)
            print ("** Debug: perl_exe" + perl_exe)
            print ("** Debug: makefile_perl: " + makefile_perl[0])

        code_perl_make = functions.system_call(perl_MakeFile_cmd)
        code_make = ""

        ##
        if (code_perl_make == 'OK'):
            ## debug messages
            if (Debug):
                print ("** Debug: perl Makefile.PL successful")

            make_bin = set_config.get_exe("make", Debug)
            print ("+ Execute make file")
            code_make = functions.system_call(make_bin)

            if (code_make == 'OK'):
                if (Debug):
                    print ("** Debug: make successful")
            else:
                print_error_message(name, package_path)
                return('n.a.')
        else:
            print_error_message(name, package_path)
            return('n.a.')
    else:
        ## perl Makefile.PL
        Build_perl = functions.retrieve_matching_files(package_path, "Build.PL", Debug)
        if (Build_perl):
            print ("+ Create Build file")
            perl_build_cmd = perl_exe + ' ' + Build_perl[0]

            ## debug messages
            if (Debug):
                print ("** Debug: chdir " + package_path)
                print ("** Debug: perl_exe" + perl_exe)
                print ("** Debug: Build_perl: " + Build_perl[0])

            code_perl_build = functions.system_call(perl_build_cmd)
            code_build = ""

            if (code_perl_build == 'OK'):
                ## debug messages
                if (Debug):
                    print ("** Debug: perl Build.PL successful")

                print ("+ Execute Build file")
                code_build = functions.system_call("./Build")

                if (code_build == 'OK'):
                    if (Debug):
                        print ("** Debug: build successful")
                else:
                    print_error_message(name, package_path)
                    return('n.a.')
            else:
                print_error_message(name, package_path)
                return('n.a.')
        else:
            print_error_message(name, package_path)
            return('n.a')


    ## copy files an finish installation
    print ("+ Include files into the PERL path")
    blib_lib = os.path.join(package_path, 'blib', 'lib')
    if os.path.isdir(blib_lib):
        return(blib_lib)
    else:
        print_error_message(name, package_path)
        return ('n.a.')
