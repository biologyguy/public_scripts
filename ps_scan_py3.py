#!/usr/bin/python3

# name: ps_scan_py3.py
# date: Nov-10-2014
# version: 1.0
# description: Python 3 PROSITE Scan client
# author: Stephen R. Bond
# email: biologyguy@gmail.com
# Repository: https://github.com/biologyguy/public_scripts
# © license: Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
# derivative work: Yes

# ======================================================================

# Derived from:
# Id: iprscan5_urllib2.py 2773 2014-04-11 11:27:27Z hpm $
# Copyright 2009-2014 EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
# ======================================================================
#
# Tested with:
#  Python 3.4.2 (Mac OS X 10.9.5)
#  !!Python 2.7.3 (Ubuntu 12.04 LTS)
#
# http://www.ebi.ac.uk/Tools/webservices/tutorials/python
# ======================================================================
# Base URL for service
baseUrl = 'http://www.ebi.ac.uk/Tools/services/rest/ps_scan'

# Load libraries
import platform, os, re, sys, time
import urllib.parse
import urllib.request
import urllib.error
from optparse import OptionParser
import xml.etree.ElementTree as eTree
from io import StringIO

# Set interval for checking status
checkInterval = 10
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

# Usage message
usage = "Usage: %prog [options...] [seqFile]"
description = """Identify protein family, domain and signal signatures in a 
protein sequence using InterProScan. For more information on InterPro and InterProScan refer to http://www.ebi.ac.uk/interpro/"""
epilog = """For further information about the InterProScan 5 (REST) web service, see http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_est."""
version = "ps_scan_py3.py 1.0 Nov-10-2014"
# Process command-line options
parser = OptionParser(usage=usage, description=description, epilog=epilog, version=version)
# Tool specific options
parser.add_option('--appl', help='signature methods to use, see --paramDetail appl')
parser.add_option('--crc', action="store_true", help='enable InterProScan Matches look-up (ignored)')
parser.add_option('--nocrc', action="store_true", help='disable InterProScan Matches look-up (ignored)')
parser.add_option('--goterms', action="store_true", help='enable inclusion of GO terms')
parser.add_option('--nogoterms', action="store_true", help='disable inclusion of GO terms')
parser.add_option('--pathways', action="store_true", help='enable inclusion of pathway terms')
parser.add_option('--nopathways', action="store_true", help='disable inclusion of pathway terms')
parser.add_option('--sequence', help='input sequence file name')
# General options
parser.add_option('--email', help='e-mail address')
parser.add_option('--title', help='job title')
parser.add_option('--outfile', help='file name for results')
parser.add_option('--outformat', help='output format for results')
parser.add_option('--async', action='store_true', help='asynchronous mode')
parser.add_option('--jobId', help='job identifier')
parser.add_option('--polljob', action="store_true", help='get job result')
parser.add_option('--status', action="store_true", help='get job status')
parser.add_option('--resultTypes', action='store_true', help='get result types')
parser.add_option('--params', action='store_true', help='list input parameters')
parser.add_option('--paramDetail', help='get details for parameter')
parser.add_option('--quiet', action='store_true', help='decrease output level')
parser.add_option('--verbose', action='store_true', help='increase output level')
parser.add_option('--baseURL', default=baseUrl, help='Base URL for service')
parser.add_option('--debugLevel', type='int', default=debugLevel, help='debug output level. Levels implemented are [1, 2, 11, 12]')

(options, args) = parser.parse_args()

if len(args) == 0:
    args = [False]

# Increase output level
if options.verbose:
    outputLevel += 1

# Decrease output level
if options.quiet:
    outputLevel -= 1

# Debug level
if options.debugLevel:
    debugLevel = options.debugLevel


# Debug print
def print_debug_message(function_name, message, level):
    if level <= debugLevel:
        print('[%s] %s' % (function_name, message), file=sys.stderr)


# User-agent for request (see RFC2616).
def get_user_agent():
    print_debug_message('get_user_agent', 'Begin', 11)
    # Agent string for urllib.request library.
    urllib_agent = 'Python-urllib/%s' % urllib.request.__version__
    client_revision = '$Revision: ???? $'
    client_version = '0'
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
        print("%s %s\n%s" % (ex.code, ex.msg, descr), file=sys.stderr)
        sys.exit()
    print_debug_message('rest_request', 'End', 11)
    return result


# Print list of parameters
def print_get_parameters():
    print_debug_message('print_get_parameters', 'Begin', 1)
    request_url = '%s/parameters' % baseUrl
    print_debug_message('print_get_parameters', 'request_url: %s' % request_url, 2)
    xml_doc = prep_xml(rest_request(request_url).decode())
    print("The following are the parameters you can set for input. Use --paramDetail <param> flag to get more details about each.")
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
    print(status.decode())
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
        print(result, file=sys.stderr)
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
        if options.outfile:
            filename = "%s.%s.%s" % (options.outfile, resultType.find("identifier").text, resultType.find("fileSuffix").text)
        else:
            filename = "%s.%s.%s" % (job_id, resultType.find("identifier").text, resultType.find("fileSuffix").text)
        # Write a result file
        if not options.outformat or options.outformat == resultType.find("identifier").text:
            # Get the result
            result = service_get_result(job_id, resultType.find("identifier").text)
            with open(filename, 'w') as fh:
                fh.write(result.decode())

            print(filename)
    print_debug_message('get_result', 'End', 1)


# Read a file
def read_file(filename):
    print_debug_message('read_file', 'Begin', 1)
    with open(filename, "r") as ifile:
        data = ifile.read().strip()
        data = re.sub("\*$", "", data)
        data = re.sub(" \t", "", data)
        sequence = re.sub(">.*\n", "", data)

    if re.search("[^A-Za-z\n]", sequence):
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
# Make sure an email address is supplied if submitting a job
elif args[0] and not options.email:
    print("Error: You must include an email address when submitting a job. E.g., $: ./iprscan5_py3.py --email YOU@EMAIL.COM my_seq_file.fasta", file=sys.stderr)
# Submit job
elif options.email and not options.jobId:
    params = {}
    if len(args) > 0:
        if os.access(args[0], os.R_OK):  # Read file into content
            params['sequence'] = read_file(args[0])
        else:  # Argument is a sequence id
            params['sequence'] = args[0]
    elif options.sequence:  # Specified via option
        if os.access(options.sequence, os.R_OK):  # Read file into content
            params['sequence'] = read_file(options.sequence)
        else:  # Argument is a sequence id
            params['sequence'] = options.sequence
    # Map flag options to boolean values.
    #if options.crc:
    #    params['crc'] = True
    #elif options.nocrc:
    #    params['crc'] = False
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
        print(new_job_id)
    else:  # Sync mode
        print(new_job_id, file=sys.stderr)
        time.sleep(5)
        get_result(new_job_id)
# Get job status
elif options.status and options.jobId:
    print_get_status(options.jobId)
# List result types for job
elif options.resultTypes and options.jobId:
    print_get_result_types(options.jobId)
# Get results for job
elif options.polljob and options.jobId:
    get_result(options.jobId)
else:
    print('Error: unrecognised argument combination', file=sys.stderr)
    parser.print_help()