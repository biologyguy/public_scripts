#!/usr/bin/env python3

"""
name: ps_scan_py3.py
date: Nov-20-2014
version: 1.0
description: Python 3 PROSITE Scan client
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/public_scripts
© license: Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
derivative work: Yes

======================================================================

Derived from:
Id: iprscan5_urllib2.py 2773 2014-04-11 11:27:27Z hpm $
Copyright 2009-2014 EMBL - European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

See:
http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
======================================================================

Tested with:
Python 3.4.2 (Mac OS X 10.9.5)
Python 3.4.0 (Linux Mint 17 LTS)

http://www.ebi.ac.uk/Tools/webservices/tutorials/python
======================================================================
"""

# Load libraries
import platform
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as eTree
from optparse import OptionParser
from io import StringIO

# Set interval for checking status
checkInterval = 10
# Output level
outputLevel = 2
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

# Usage message
usage = "Usage: %prog [options...] [seqFile]"
description = """Identify protein family, domain and signal signatures in a protein sequence using InterProScan.
For more information on InterPro and InterProScan refer to http://www.ebi.ac.uk/interpro/"""
epilog = """For further information about the PROSITE Scan (REST) web service, see
http://www.ebi.ac.uk/Tools/webservices/services/pfa/ps_scan_rest"""
version = "ps_scan_py3.py 1.0 Nov-20-2014"
# Process command-line options
parser = OptionParser(usage=usage, description=description, epilog=epilog, version=version)
# Tool specific options
parser.add_option('--appl', help='signature methods to use, see --paramDetail appl')
parser.add_option('--goterms', action="store_true", help='enable inclusion of GO terms')
parser.add_option('--nogoterms', action="store_true", help='disable inclusion of GO terms')
parser.add_option('--pathways', action="store_true", help='enable inclusion of pathway terms')
parser.add_option('--nopathways', action="store_true", help='disable inclusion of pathway terms')
parser.add_option('--sequence', action="store", help='Input sequence explicitly on command line')
# General options
parser.add_option('--email', help='e-mail address')
parser.add_option('--title', help='job title')
parser.add_option('--outfile', help='file name for results')
parser.add_option('--outformat', help='output format for results')
parser.add_option('--async', action='store_true', help='asynchronous mode')
parser.add_option('--jobId', action="store", help='job identifier')
parser.add_option('--polljob', action="store_true", help='get job result')
parser.add_option('--status', action="store_true", help='get job status')
parser.add_option('--resultTypes', action='store_true', help='get result types')
parser.add_option('--params', action='store_true', help='list input parameters')
parser.add_option('--paramDetail', help='get details for parameter')
parser.add_option('--outputLevel', type=int,
                  help='Explicilty set the output verbosity. 0 == quiet, 3 == verbose, 1 and 2 are intermediate.')
parser.add_option('--quiet', action='store_true', help='decrease output level')
parser.add_option('--verbose', action='store_true', help='increase output level')
parser.add_option('--service', choices=['prosite_scan', 'interpro'], default='prosite_scan',
                  help='Which EMBL-EBI REST service do you want?')
parser.add_option('--debugLevel', type='int', default=debugLevel,
                  help='debug output level. Levels implemented are [1, 2, 11, 12]')

(options, args) = parser.parse_args()

# Base URL for service
if options.service == "prosite_scan":
    baseUrl = 'http://www.ebi.ac.uk/Tools/services/rest/ps_scan'
elif options.service == "interpro":
    baseUrl = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5'


if len(args) == 0:
    args = [False]

# Set output level (note order of precedence)
if options.outputLevel:
    outputLevel = options.outputLevel

if options.verbose:
    outputLevel = 3

if options.quiet:
    outputLevel = 0

# Debug level
if options.debugLevel:
    debugLevel = options.debugLevel


# Debug print
def print_debug_message(function_name, message, level):
    if level <= debugLevel:
        print('[%s] %s' % (function_name, message), file=sys.stderr)


def print_stdout(message, level, line_break=True):
    if level <= outputLevel:
        if line_break:
            print(message)
        else:
            print(message, end="")


# User-agent for request (see RFC2616).
def get_user_agent():
    print_debug_message('get_user_agent', 'Begin', 11)
    # Agent string for urllib.request library.
    urllib_agent = 'Python-urllib/%s' % urllib.request.__version__
    client_revision = '$Revision: ???? $'
    client_version = '1.0'
    if len(client_revision) > 11:
        client_version = client_revision[11:-2]
    # Prepend client specific agent string.
    user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        client_version, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        urllib_agent
    )
    print_debug_message('get_user_agent', 'user_agent: %s' % user_agent, 12)
    print_debug_message('get_user_agent', 'End', 11)
    return user_agent


# Wrapper for a REST (HTTP GET) request
def rest_request(url):
    print_debug_message('rest_request', 'Begin', 11)
    print_debug_message('rest_request', 'url: %s' % url, 11)
    # Errors are indicated by HTTP status codes.
    try:
        # Set the User-agent.
        user_agent = get_user_agent()
        http_headers = {'User-Agent': user_agent}
        req = urllib.request.Request(url, None, http_headers)
        # Make the request (HTTP GET).
        req_h = urllib.request.urlopen(req)
        result = req_h.read()
        req_h.close()
    # Errors are indicated by HTTP status codes.
    except urllib.error.HTTPError as ex:
        # Trap exception and output the document to get error message.
        error = prep_xml(ex.file.read().decode())
        descr = error.find("description").text
        print("%s %s\n%s" % (ex.code, ex.msg, descr))
        sys.exit()
    print_debug_message('rest_request', 'End', 11)
    return result


# Print list of parameters
def print_get_parameters():
    print_debug_message('print_get_parameters', 'Begin', 1)
    request_url = '%s/parameters' % baseUrl
    print_debug_message('print_get_parameters', 'request_url: %s' % request_url, 2)
    xml_doc = prep_xml(rest_request(request_url).decode())
    print("The following are the parameters you can set for input. Use --paramDetail <param> flag to get "
          "more details about each.")
    for next_id in xml_doc.findall("id"):
        print("\t%s" % next_id.text)
    print_debug_message('print_get_parameters', 'End', 1)    


# Print description of a parameter
def print_get_parameter_details(param_name):
    print_debug_message('print_get_parameter_details', 'Begin', 1)
    request_url = '%s/parameterdetails/%s' % (baseUrl, param_name)
    xml_tree = prep_xml(rest_request(request_url).decode())
    param_name = xml_tree.find("name").text
    param_type = xml_tree.find("type").text
    param_desc = xml_tree.find("description").text
    print("%s\t(type: %s)" % (param_name, param_type))
    print("%s\n" % param_desc)
    values = xml_tree.find("values")
    if values:
        print("Valid values accepted:")
        for value in values.findall("value"):
            label = value.find("label").text
            val = value.find("value").text
            default = value.find("defaultValue").text
            output = ""
            if val:
                output += val.ljust(15)
                output += "---\t"
                if label:
                    output += "%s" % label
                if default == "true":
                    output += "  (default)"
            properties = value.find("properties")
            if properties:
                for prop in properties:
                    key = prop.find("key").text
                    prop_val = prop.find("value").text
                    if key:
                        output += "\n\t%s:" % key
                    if prop_val:
                        output += "\t%s" % prop_val
            print(output)
    print_debug_message('print_get_parameter_details', 'End', 1)


# Submit job
def service_run(email, title, run_params):
    print_debug_message('service_run', 'Begin', 1)
    # Insert e-mail and title into params
    run_params['email'] = email
    if title:
        run_params['title'] = title
    request_url = '%s/run/' % baseUrl
    print_debug_message('service_run', 'request_url: %s' % request_url, 2)
    # Signature methods requires special handling (list)
    appl_data = ''
    if 'appl' in run_params:
        # So extract from params
        appl_list = run_params['appl']
        del run_params['appl']
        # Build the method data options
        for appl in appl_list:
            appl_data += '&appl=%s' % appl
    # Get the data for the other options
    request_data = urllib.parse.urlencode(run_params)
    # Concatenate the two parts.
    request_data += appl_data
    print_debug_message('service_run', 'requestData: %s' % request_data, 2)
    # encode requestData for transmission
    request_data = request_data.encode()

    # Errors are indicated by HTTP status codes.
    try:
        # Set the HTTP User-agent.
        user_agent = get_user_agent()
        http_headers = {'User-Agent': user_agent}
        req = urllib.request.Request(request_url, None, http_headers)
        # Make the submission (HTTP POST).

        req_h = urllib.request.urlopen(req, request_data)
        job_id = req_h.read().decode("utf-8")
        req_h.close()

    except urllib.error.HTTPError as ex:
        # Trap exception and output the document to get error message.
        print(ex.read(), file=sys.stderr)
        raise
    print_debug_message('service_run', 'job_id: %s' % job_id, 2)
    print_debug_message('service_run', 'End', 1)
    return job_id


# Get job status
def service_get_status(job_id):
    print_debug_message('service_get_status', 'Begin', 1)
    print_debug_message('service_get_status', 'job_id: %s' % job_id, 2)
    request_url = '%s/status/%s' % (baseUrl, job_id)
    print_debug_message('service_get_status', 'request_url: %s' % request_url, 2)
    status = rest_request(request_url)
    print_debug_message('service_get_status', 'status: %s' % status.decode(), 2)
    print_debug_message('service_get_status', 'End', 1)
    return status


# Print the status of a job
def print_get_status(job_id):
    print_debug_message('print_get_status', 'Begin', 1)
    status = service_get_status(job_id)
    print_stdout(status.decode(), 2)
    print_debug_message('print_get_status', 'End', 1)
    

def prep_xml(xml):
    file_like = StringIO(xml)
    xml_tree = eTree.parse(file_like)
    return xml_tree.getroot()


# Get available result types for job
def service_get_result_types(job_id):
    print_debug_message('service_get_result_types', 'Begin', 1)
    print_debug_message('service_get_result_types', 'job_id: %s' % job_id, 2)
    request_url = '%s/resulttypes/%s' % (baseUrl, job_id)
    print_debug_message('service_get_result_types', 'request_url: %s' % request_url, 2)
    xml_doc = prep_xml(rest_request(request_url).decode())
    output = []
    for child in xml_doc:
        output.append(child)
    print_debug_message('service_get_result_types', 'End', 1)
    return output


# Print list of available result types for a job.
def print_get_result_types(job_id):
    print_debug_message('print_get_result_types', 'Begin', 1)
    result_type_list = service_get_result_types(job_id)
    for resultType in result_type_list:
        identifier = resultType.find("identifier")
        print(identifier.text)
        label = resultType.find("label").text
        descr = resultType.find('description').text
        media_type = resultType.find('mediaType').text
        file_suffix = resultType.find('fileSuffix').text
        if label:
            print("\tLabel:\t\t", label)
        if descr:
            print("\tDescription:\t", descr)
        if media_type:
            print("\tMedia type:\t", media_type)
        if file_suffix:
            print("\tFile suffix:\t", file_suffix)
    print_debug_message('print_get_result_types', 'End', 1)


# Get result
def service_get_result(job_id, result_type):
    print_debug_message('service_get_result', 'Begin', 1)
    print_debug_message('service_get_result', 'job_id: %s' % job_id, 2)
    print_debug_message('service_get_result', 'type: %s' % result_type, 2)
    request_url = '%s/result/%s/%s' % (baseUrl, job_id, result_type)
    result = rest_request(request_url)
    print_debug_message('service_get_result', 'End', 1)
    return result


# Client-side poll
def client_poll(job_id):
    print_debug_message('client_poll', 'Begin', 1)
    result = 'PENDING'
    while result == 'RUNNING' or result == 'PENDING':
        result = service_get_status(job_id).decode("utf-8")
        print_stdout(result, 2)
        if result == 'RUNNING' or result == 'PENDING':
            time.sleep(checkInterval)
    print_debug_message('client_poll', 'End', 1)


# Get result for a job_id
def get_result(job_id):
    print_debug_message('get_result', 'Begin', 1)
    print_debug_message('get_result', 'job_id: %s' % job_id, 1)
    # Check status and wait if necessary
    client_poll(job_id)
    # Get available result types
    result_types = service_get_result_types(job_id)
    for resultType in result_types:
        # Derive the filename for the result
        identifier = resultType.find("identifier").text
        file_suffix = resultType.find("fileSuffix").text
        if options.outfile:
            filename = "%s.%s.%s" % (options.outfile, identifier, file_suffix)
        else:
            filename = "%s.%s.%s" % (job_id, identifier, file_suffix)
        # Write a result file
        if not options.outformat or options.outformat == identifier:
            # Get the result
            result = service_get_result(job_id, identifier)
            with open(filename, 'wb') as fh:
                fh.write(result)
            print_stdout(filename, 3)
    print_debug_message('get_result', 'End', 1)


# Read a file
def read_file(filename):
    print_debug_message('read_file', 'Begin', 1)
    with open(filename, "r") as ifile:
        data = ifile.read().strip()
        data = re.sub("\*$", "", data)
        data = re.sub(" \t", "", data)

    seq = re.sub(">.*\n", "", data)  # The server will handle FASTA headers, but no funny characters in sequence
    if re.search("[^A-Za-z\n]", seq):
        sys.exit("Error: Invalid characters found in the sequence provided.")

    print_debug_message('read_file', 'End', 1)
    return data

# No options... print help.
if numOpts < 2:
    parser.print_help()

# List parameters
elif options.params:
    print_get_parameters()

# Get parameter details
elif options.paramDetail:
    print_get_parameter_details(options.paramDetail)


# Get job status
elif options.status:
    if not options.jobId:
        sys.exit("Error: You must include --jobId to retrieve status.")
    else:
        print_get_status(options.jobId)

# List result types for job
elif options.resultTypes:
    if not options.jobId:
        sys.exit("Error: You must include --jobId to retrieve result types.")
    else:
        print_get_result_types(options.jobId)

# Get results for job
elif options.polljob:
    if not options.jobId:
        sys.exit("Error: You must include --jobId to poll a job.")
    else:
        get_result(options.jobId)

# Submit job
elif args[0] or options.sequence:
    # Make sure an email address is supplied if submitting a job
    if not options.email:
        sys.exit("Error: You must include an email address when submitting a job. E.g., $: ./iprscan5_py3.py --email "
                 "YOU@EMAIL.COM my_seq_file.fasta")

    params = {}
    if args[0]:
        if os.access(args[0], os.R_OK):  # Read file into content
            params['sequence'] = read_file(args[0])
        else:  # Argument is a sequence id
            params['sequence'] = args[0]
    elif options.sequence:  # Passing in the actual sequence on command line
        sequence = options.sequence.strip()
        sequence = re.sub("\*$", "", sequence)
        sequence = re.sub(" \t", "", sequence)
        if re.search("[^A-Za-z\n]", re.sub(">.*\n", "", sequence)):
            sys.exit("Error: Invalid characters found in the sequence provided.")
        params['sequence'] = sequence

    # Map flag options to boolean values.
    if options.goterms:
        params['goterms'] = True
    elif options.nogoterms:
        params['goterms'] = False
    if options.pathways:
        params['pathways'] = True
    elif options.nopathways:
        params['pathways'] = False
    # Add the other options (if defined)
    if options.appl:
        params['appl'] = re.split('[ \t\n,;]+', options.appl)
    
    # Submit the job
    new_job_id = service_run(options.email, options.title, params)
    if options.async:  # Async mode
        print_stdout("Project ID: ", 2, line_break=False)
        print_stdout(new_job_id, 1)
    else:  # Sync mode
        print_stdout("Project ID: ", 2, line_break=False)
        print_stdout(new_job_id, 1)
        time.sleep(5)
        get_result(new_job_id)

else:
    print('Error: unrecognised argument combination', file=sys.stderr)
    parser.print_help()
