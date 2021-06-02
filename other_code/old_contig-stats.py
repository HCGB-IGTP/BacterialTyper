def old_contig_stats(assembly_file, csv_arg):
    """Generate assembly statistics
    
    Calls additional perl script to generate contig statistics. Steps:
    
    - Retrieve information of additional perl script location (using :func:`BacterialTyper.other_tools.tools.perl_scripts`)
    
    - Create system call (using :func:`BacterialTyper.scripts.functions.system_call`) and return output statistic file created.
    
    :param assembly_file: Absolute path to assembly fasta file.
    :type assembly_file: string
    :param csv_arg: Comma separated values for splitting sets. Default: 1000,10000
    :type csv_arg: string
    :return: Text file containing statistics for later analysis.
    :rtype: string
    
    The perl script contig_stats.pl has a mandatory argument which is a single fasta file and default splitting sets are: 1000, 10000. User can provide new parts using a csv argument for the script:
    
    .. code-block:: sh

        perl BacterialTyper/other_tools/perl/contig_stats.pl fasta_file 1000,10000
        

    .. seealso:: This function depends on other BacterialTyper functions called:
    
        - :func:`BacterialTyper.scripts.functions.system_call`
        
        - :func:`BacterialTyper.other_tools.tools.perl_scripts`
    """
    contig_stats_script = tools.perl_scripts('contig_stats')

    
    file_out = assembly_file + '_stats.txt'
    cmd_stats = 'perl %s %s %s > %s' %(contig_stats_script, assembly_file, csv_arg, file_out) ## [TODO] Generate this code in python
    code_chr = functions.system_call(cmd_stats)
    return (file_out)
