#!/usr/bin/env python
"""
Use this script to serve files from the login nodes.

We can use this script to set up a temporary web server on the examples directory,
and view any generated images directly in a browser. It serves files under the
current directory where the script is executed. Use with care.
"""

import sys
import socket
from SimpleHTTPServer import SimpleHTTPRequestHandler as SimpleHandler
from SocketServer import TCPServer

# usage
if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
    print "Usage: {0} [PORT]".format(sys.argv[0])
    sys.exit(1)

# set the port number
PORT = 8000
if len(sys.argv) > 1:
    PORT = int(sys.argv[1])

# use a default mime type of 'text/plain', rather than 'application/octet-stream'
extensions = ['', '.sh', '.rst']
for ext in extensions:
    SimpleHandler.extensions_map[ext] = 'text/plain'

# create a tcp server on the given port
httpd = TCPServer(("", PORT), SimpleHandler)

# run the server
print "On http://{0}.hoffman2.idre.ucla.edu:{1}".format(socket.gethostname(), PORT)
print "Serving HTTP on 0.0.0.0 port {0}".format(PORT)
httpd.serve_forever()

# EOF
