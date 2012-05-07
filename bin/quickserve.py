#!/usr/bin/env python
"""
Use this script to serve files from the login nodes of the cluster.

We can use this script to set up a temporary web server on the examples directory,
and view any generated images directly in a browser. It serves files under the
current directory where the script is executed. Use with care.
"""

import sys
import socket
import base64
from SimpleHTTPServer import SimpleHTTPRequestHandler as SimpleHandler
from SocketServer import TCPServer

# command line arguments
from optparse import OptionParser
parser = OptionParser(usage="Usage: %prog [-p PORT] [-a LOGIN:PASSWORD]")
parser.add_option('-p', '--port', help='port number')
parser.add_option('-a', '--auth', help='authentication of the form LOGIN:PASSWORD')
(options, args) = parser.parse_args()

# set the port number
PORT = 8000
if options.port:
    PORT = int(options.port)

# set the auth information (colon-separated login and password)
AUTH = 'sloc:test'
if options.auth:
    AUTH = options.auth
LOGIN = AUTH.split(':')[0]
PASS = AUTH.split(':')[1]
BASIC_AUTH = 'Basic ' + base64.b64encode(AUTH)

# extend SimpleHandler with basic access authentication
class Handler(SimpleHandler):
    """
    Add basic access authentication to the SimpleHTTPRequestHandler.
    http://en.wikipedia.org/wiki/Basic_access_authentication
    """
    def do_AUTHHEAD(self):
        self.send_response(401)
        self.send_header('WWW-Authenticate', 'Basic realm="quickserve.py"')
        self.send_header('Content-type', 'text/html')
        self.end_headers()
    def do_GET(self):
        auth = self.headers.getheader('Authorization')
        if auth is None:
            self.do_AUTHHEAD()
            self.wfile.write('no auth header received')
            return
        elif auth == BASIC_AUTH:
            SimpleHandler.do_GET(self)
            return
        else:
            self.do_AUTHHEAD()
            self.wfile.write('not authenticated!')
            return

# use a default mime type of 'text/plain', rather than 'application/octet-stream'
extensions = ['', '.sh', '.rst']
for ext in extensions:
    Handler.extensions_map[ext] = 'text/plain'

# create a tcp server on the given port
httpd = TCPServer(("", PORT), Handler)

# run the server
print "On http://{0}.hoffman2.idre.ucla.edu:{1}".format(socket.gethostname(), PORT)
print "Using login '{0}' and password '{1}'".format(LOGIN, PASS)
print "Serving HTTP on 0.0.0.0 port {0}".format(PORT)
httpd.serve_forever()

# EOF
