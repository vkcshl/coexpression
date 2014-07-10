#!/usr/bin/python
import argparse
import sys
import os
import time
import traceback
import sys
import ctypes
import subprocess
import os
from optparse import OptionParser


desc1 = '''
NAME
      coex-filter-genes -- select differentially expressed genes

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Check status of probabilistic annotation jobs submitted by the user.  For
      each job, information about the job is displayed.  A job that has completed
      is then deleted from the system.

      The --jobID optional argument is the identifier of a specific job to
      check.

      The ujs-url optional argument specifies an alternate URL for the user and
      job state service.
'''

desc3 = '''
EXAMPLES
      Filter genes with ANOVA
      > coex-filter-genes
      
      Filter genes with LOR
      > coex-filter-genes 

SEE ALSO
      coex_filter

AUTHORS
      
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='coex-filter-genes', epilog=desc3)
    parser.add_argument('-u', '--ws_url', help='Workspace url', action='store', dest='ws_url', default='https://kbase.us/services/ws')
    parser.add_argument('-w', '--ws_id', help='Workspace id', action='store', dest='ws_id', default=None, required=True)
    parser.add_argument('-i', '--in_id', help='Input Series object id', action='store', dest='inobj_id', default=None, required=True)
    parser.add_argument('-o', '--out_id', help='Output Series object id', action='store', dest='outobj_id', default=None, required=True)
    parser.add_argument('-n', '--num_genes', help='The number of genes to be selected', action='store', dest='num_genes', default=None)
    parser.add_argument('-p', '--p_value', help='The p-value cut-off', action='store', dest='p_value', default=None)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    if(args.p_value is None and args.num_genes is None) :
        print "Either p_value or num_genes has to be specified";
        exit(1);


    exit(0);
    
    # Get the list of jobs for the user.
    if 'KB_AUTH_USER_ID' in os.environ:
        userID = os.environ.get('KB_AUTH_USER_ID')
    else:
        auth = _read_inifile()
        userID = auth['user_id']
    ujsClient = UserAndJobState(args.ujsURL)
    try:
        jobList = ujsClient.list_jobs([ userID ], 'RCE')
    except JobStateServerError as e:
        print e.message
        exit(1)
    
    # See if the user has any jobs in the list.
    if len(jobList) == 0:
        print 'There are no jobs for you.'
        exit(1)

    # Print info about the specific job if requested.
    if args.jobID is not None:
        for job in jobList:
            info = job_info_dict(job)
            if args.jobID == info['id']:
                print_job(info)
                exit(0)
        print "Job '%s' was not found." %(args.jobID)
        exit(1)
        
    # Print all of the jobs in the list.
    for job in jobList:
        info = job_info_dict(job)
        if info['description'].split()[0] == 'pa-annotate': # Only show pa-annotate jobs
            print_job(info)
            
    exit(0)
